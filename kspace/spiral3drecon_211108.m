function [tmpk, im_all] = spiral3drecon_211108(myPfile,doRots, Ndel,  doNUFFT)

%function [tmpk im_all] = spiral3drecon_201119(myPfile, doRots, Ndel, doNUFFT)
%
% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot - USES rotation matrices sotred in text files.
% 3 - reconstructs the image by Gridding and by Jeff's NUFFT
%
% Inputs:
%   myPfile: P file name
%   doRots: is it rotating (1), or stack of spirals(0)
%   Ndel : number of samples to shift the FID
%   doNUFFT:  use the NUFFT reconstruction from the MIRT toolbox

dt = 4e-6; %seconds
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

% configuration for reconstruction:
SHOWPIX = 3 ;     % DEBUG and show pix of the different steps
scaleFactor = 1;
doT2filter = 0;
T2 = 0.150;
Nramp = 20;         % Number of ramp points to throw out
doPhaseCyle = 0;   % CPMG phase cycling
doPhaseDetrend =1;%   linear phase correction for each echo
doYrot = 0;
doGoldenAngle = 1; %   doGoldenAngle: rotations done by Golden Angle (137.5) vs. 180/n.slices?
useKrot = 0;
Ncenter = 0;
zoomfactor = 1.5;

if nargin==0
    pfiles = dir('P*.7');
    myPfile = pfiles(1).name;
    doNUFFT = 0;
    doRots = 1;
    Ndel = -3;
    Ndel2 = Ndel;
    
end
Ndel2 = Ndel ;


if doRots == 1
    useKrot = 1; %   useKrot:  read in krots.txt and use it to calculate the rotations from interleaves v. slice encodes
else
    useKrot = 0;
    %doPhaseDetrend = 0;
end

fprintf('\nBegin recon of 3D spiral trajectory for %s, delay=%d ...', myPfile, Ndel)
%run('~hernan/matlab/irt/setup')

% load the raw data from P file
% Nframes x Nslices*ndat x Ncoils
[raw scaninfo ] = read_raw_3d(myPfile ,0 );

% parse out the header information
nslices = scaninfo.nslices;
nechoes = scaninfo.npr;
nframes = scaninfo.nphases;
ndat = scaninfo.ndat;
ncoils = scaninfo.ncoils;
nfids = nechoes*nframes*nslices;
fov = scaninfo.opfov;
MYDIM = scaninfo.opxres;
nleaves = scaninfo.npr;
slthick = scaninfo.slthick/10;
fov3 = [fov fov nslices*slthick];
dim3 = [MYDIM MYDIM nslices];
MaxFrame = nframes;
%%%% for debugging:
% MaxFrame = 1;


if doRots
    fov3 = [fov fov fov];
    dim3 = [MYDIM MYDIM MYDIM];
end
dim3 = round(dim3/zoomfactor);
dim3 = dim3+1;  % the center of k-space is important, so make sure the dimensions allow it

% make some headers for output
h=define_nii_hdr();
h.dim = [4 dim3 nframes 0 0 0];
h.pixdim = [1 fov3./dim3 1 0 0 0];
h.datatype = 4;
h.bitsperpixel = 16;
ah = nii2avw_hdr(h);

if SHOWPIX >2
    figure(1)
    subplot(221)
    plot(squeeze((abs(raw(1,:,:)))));
    title('Echo train: Frame 1, all coils');
    %fpos = get(gcf,'Position');
    %set(gcf,'Position', fpos + [0 600 0 0]);
    drawnow
end

echo = squeeze((abs(raw(1,:,:))));


% First get the k-space trajectory configured:

fprintf('\nLoading grads...');
g = load('grad.txt');

glen = size(g,1);
aqwindow = dt*glen;
dt2 = aqwindow/ndat;
cnt = 1;


% Calculate the k-space trajectory calculated from the stored 
% gradient waveformss, rather than assuming an ideal trajectory
k2 = gamma*cumsum(g)*dt;  % cm-1 units
k2 =  k2';

% calculate the k-space boundaries of the readout:
Kxymax = MYDIM/fov/2;  % Kmax
Kzmax = Kxymax;

% must resample the trajectory to match the data acquisition rate
ks = zeros(3,ndat);
fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
for n=1:3
    tmp = k2(n,:);
    ks(n,:) = interp1(linspace(0,1,length(tmp)) , tmp, linspace(0,1,ndat));
end

% throw out ramp points if specified
ks = ks(:, Nramp+1:end-Nramp) ;

% 11.20.21 : throw out center points
ks(:, end/2-Ncenter/2:end/2+Ncenter/2) = [];

% 11.21.21 : add a different weight to the spiral out ?
Nsamp = size(ks,2);
weights = ones(1,Nsamp);
% weights(end/2+1:end) = i;
% I can add a delay here as a linear phase?
%weights = exp(-i * Ndel * 2 * pi  * linspace(0,1, Nsamp) /Nsamp);

if doRots  % SERIOS readout: rotating spiral pancakes about x axis
        
    ks_all = [];
    echolist = [];
    idx=0;
    
    
    phi_increment = pi/nleaves;
    ks_plate = -ks;
    
    for leafn = 0:nleaves-1
        for s = 0:nslices-1
            glen = size(ks_plate,1);
            idx = leafn*nslices + s +1 ;
            if useKrot
                % instead of calculating and  building our own trajectory,
                % use whatever rotating matrices are written out by the scanner
                % in the text file "krotations.txt"
                %
                %     The list of rotations is written in this order
                %     (same as acquisition):
                %     leaf 1:  slice 1 , slice 2 ... leaf 2:  slice 1 ,slice 2...
                %     NOTE: However, the data are in arranged in this order:
                %     slice 1: leaf 1, leaf 2 ... slice 2:  ...etc
                %
                % The first four columns in the file indicate:
                % slice #, leaf #, x-rotations angle, z-rotation angle
                %
                % the last 9 columns are the actual rotation matrix (redundant)
                
                krot = load('krotations.txt');
                angle_list = krot(:,3:4);
                
                % not used after all:
                krot = krot(:, 5:end);
                krot = krot/32767;

                row = krot(idx,:);
                rot = reshape(row, 3,3)'; 
                tmpks = rot*ks_plate;
                tmpks = tmpks;
                
                % use the angles from the file
                % and generate rot matrices from those
                % rotation angles are in the first two columns of the file
%                 xi = angle_list(idx,1);
%                 phi = angle_list(idx,2);
                
            else
                % calculate rotation angles on the fly,
                % assuming GoldenAngle rotations about x-axis and (slices)
                % evenly spaced rotations about z-axis (leaves)
                
                % xi : angle for Rx and Ry rotations
                xi = deg2rad(137.5)* (s + leafn/nleaves);  % slice rotations (xi)
                % phi  : angle for Rz rotations
                phi = phi_increment*leafn ;      % leave rotations  (phi)
                
                
                % 8.6.2020: the rotation matrix is doing (-Rx)*Rz
                % xi needs to be negative in order for rotations to agree
                % with angles
                sxi = sin(xi); cxi = cos(xi);
                sphi = sin(phi); cphi = cos(phi);
                
                % Note that this gives the right orientation at the end, 
                % although it's
                % different from krots.
                Rx = [1  0   0 ;
                    0  cxi sxi ;
                    0 -sxi cxi];
                
                Ry = [cxi  0   sxi ;
                      0    1   0 ;
                     -sxi  0 cxi];

                Rz = [cphi -sphi 0;
                     sphi cphi 0;
                     0   0   1];
                
                if ~doYrot
                    Ry = eye(3);
                end
                
                tmpks = Rz * ks_plate;
                tmpks = Rx * tmpks;
                tmpks = Ry * tmpks;

                rotmat = (Ry*Rx*Rz)';
                
                echolist = [echolist; s, leafn,  xi, phi,  32767*rotmat(:)'];
                
                
            end
            if SHOWPIX >3
                
                figure(1)
                subplot(2,2,2)
                plot3(tmpks(1,:), tmpks(2,:), tmpks(3,:));
                axis(2*Kxymax*[-1 1 -1 1 -1 1]);
                hold off
                drawnow
                
                
            end
            
            ks_all = [ks_all; tmpks'];
            
            
    
        end
    end
    
    ks = ks_all;
    
    if SHOWPIX >2
        figure(1)
        subplot(2,2,2)
        plot3(ks(:,1), ks(:,2), ks(:,3));
        axis(2*Kxymax*[-1 1 -1 1 -1 1]);
        hold off
        
       
    end
    
    
    
    
    
else  % stack of spirals
    
    % calculate the kz components for each echo
    Kzmax = 1/(slthick)/2;
   
    % the kz steps order and corresponding fraction of Gz max are
    % specified in a file
    kzsteps = load('kzsteps.txt');
    kzfraction = -kzsteps(:,3);      % faction of kz space for each platter
    kzorder = kzsteps(:,2) ;      % order of the platters during readout  (center -> out)
%     
%     for n=1:numel(kzfraction)
%         kzfraction(kzorder(n)) = kzsteps(n,3);      % faction of kz space for each platter
%     end
%  
%    kzfraction = -kzfraction;
%      
    % calculate kx, ky kz locations for each platter:
    % kx and ky are always the same, but kz needs to 
    % be scaled for every platter
    
    % copy the xy trajectory to all platters
    ks_plate = ks;
    ks_all = [];
    
    % now re-scale:
    fin = 0;
    % if  more than one interleave per platter,
    % we need to rotate about the z-axis for each
    % interleaf.
                    
    for m=1:nleaves

        for n=1:nslices

            phi = (m-1)*pi/nleaves;
            sphi = sin(phi); 
            cphi = cos(phi);
            
            Rz = [cphi -sphi 0;
                sphi cphi 0;
                0 0 1];
            
            Ktmp = Rz * ks_plate; 
            Ktmp(3,:) = Ktmp(3,:) * kzfraction(n);

            % the kz is not rotation, but increasing by kzfraction
            ks_all = [ks_all; Ktmp']; 

            if SHOWPIX >2
                figure(1)
                subplot(2,2,2)
                plot3(Ktmp(1,:), Ktmp(2,:), Ktmp(3,:));
                axis(2*Kxymax*[-1 1 -1 1 -1 1]);
                hold off
                pause(0.1)
                drawnow
            end
            
        end
    end
    
    % rejoin the coordinates of the trajectory into a single matrix
    ks = ks_all;
       
    % 2D test
    % kzs(beg:fin) = 0;
    
end
weights = repmat(weights,1,nslices*nechoes);
if SHOWPIX >3
 
    figure(6)
    subplot(311), plot(diff(ks(:,1)))
    subplot(312), plot(diff(ks(:,2)))
    subplot(313), plot(diff(ks(:,3)))
    drawnow
    
end

fprintf('\nk-space trajectory ready...\n' );
  
% k-space Trajectory ready .
% now process the data and recon the images

% the dimensions of raw are: Nframes x Nslices*ndat x Ncoils
imseries = zeros(dim3(1)*dim3(2)*dim3(3),nframes);

% we need to generate a cartesian interpolation
% grid for the data

[kx , ky, kz ]= meshgrid( ...
    linspace(-Kxymax,Kxymax,dim3(1)), ...
    linspace(-Kxymax,Kxymax,dim3(2)), ...
    linspace(-Kzmax,Kzmax,dim3(3)));


if SHOWPIX > 1
    MaxFrame = 1;
end

raw = raw * scaleFactor;


fprintf('\nMaking Density compensation function...');
% distance to center of kspace
dist= ((ks(:,1).^2 + ks(:,2).^2 + ks(:,3).^2 )) .^0.5;
kdens = ones(size(dist));
if doNUFFT
    
    %w = (tanh(linspace(-pi, pi, length(dist))) +1 )/2;
    %kdens = dist .* w';
    
    % density compensaton function :number of samples in a region within
    % 2*deltaK from each sample.
    % kdens = (sampdensity_dist(ks, Kzmax/(MYDIM/2))).^(-1);
    % (we get a sharper image  if we use bigger exponent)
    %kdens(isnan(kdens)) = 0;
    %kdens = 1e-4 + kdens/max(kdens);
    
    % OR ... use the average distance to the nearest 6 neighbor smaples
    kdens = sampdensity(ks, 6);
    
    % kdens = 1e-3 + kdens/max(kdens);
    % smoother image if we use larger DC term
    save dcf.mat kdens dist
end

if SHOWPIX==5
    figure(4)
    plot(kdens(1:ndat-2*Nramp));
    title('K-density compensation');
    
end


for f=1 :MaxFrame
    
    im_all = 0;
    coilbuf = zeros(dim3(1), dim3(2), dim3(3), ncoils);
    
   parfor c = 1:ncoils
%   for c = 1:ncoils
        fprintf('\nReconning frame %d , coil %d, delay %d...\n', f, c, Ndel);
        
        % reshape, filter and adjust the data:
        
        tmpk = raw(f,:,c);   % new dims:  Nslices*ndat*nleaves x 1
        %
        tmpk = reshape(tmpk, ndat*nleaves, nslices);
        %tmpk = reshape(tmpk, ndat*nslices, nleaves);
        
         
        b = reshape(tmpk,ndat,nleaves,nslices);
        
    
        % accouting for CPMG phase cycling during acquisition:
        if doPhaseCyle
            
            if doRots % SERIOS readout
                for ss = 1:2:nslices
                    b(:,:,ss) =  exp(-i*pi)*(b(:,:,ss))  ;
                end
                
            else  % stack of spirals:  different order
                
                for ss = 1:nslices/2
                    b(:,:,ss) =  exp(-i*pi)*(b(:,:,ss))  ;
                end
            end
            %{            
            % 11/19/21   reverse every other echo instead of above?
            b(:,1:2:end,:) = b(end:-1:1, 1:2:end, :)  ;
            %}
        end
        
        % Remove linear phase trend in the data
        if doPhaseDetrend
            for shot = 1:nechoes
                for ss = 1:nslices
                    ph_ref = angle(b(round(end/2), 1, ss));
                    b(:,shot,ss) = phaseDetrend(b(:,shot,ss));
                   
                    if doRots == 0
                         b(:,shot,ss) = phaseDetrend2(b(:,shot,ss), ph_ref);
                    end
                    % make the phase at the center of all interleaves match
                    %ph_correction = angle(b(end/2, shot, ss)) - angle(b(end/2, 1, 1)) ; 
                    %b(:,shot, ss) = b(:,shot, ss) .* exp(-i * ph_correction);
                end
            end
            
        end
        
        % add a T2 filter
        if doT2filter
            if ~doRots
                % the reference is the center of the first echo (averaged
                % over leaves)
                % rescale other echoes so that the center of k-space has
                % the same amplitude as the first echo
                kcenter0 = mean(abs(b(end/2-2:end/2+2, 1,1)));
                
                for ss = 1:nslices
                    for ll=1:nleaves
                        kcenter = mean(abs(b(end/2-2:end/2+2, ll,ss)));
                        b(:,ll,ss) = b(:,ll,ss) * kcenter0 / kcenter;
                    end
                end
                
            else
                % the reference is a few points at the center ot the first echo
                % rescale other echoes so that the center of k-space has
                % the same amplitude as the first echo
                kcenter0 = mean(abs(b(end/2-2:end/2+2, 1,1)));  
                for ss = 1:nslices
                    for ll=1:nleaves
                        kcenter = mean(abs(b(end/2-2:end/2+2, ll,ss)));
                        b(:,ll,ss) = b(:,ll,ss) * kcenter0 / kcenter;
                    end
                end
            end
        end
        
        %b  = circshift(b, [Ndel 0 0]);
        % include delays on each half of the spiral separately

        bin  = circshift(b, [Ndel 0 0]);
        bout = circshift(b, [Ndel2 0 0]);

        bin = bin(1:end/2, :,:);
        bout = bout(end/2+1:end, :,:);

        % fix the ends
        for cnt=1:size(b,2)
            bin(1:abs(Ndel), cnt) = mean(bin(10:15,cnt),1);
            bout(end-abs(Ndel):end, cnt) = mean(bout(end-10:end-15,cnt),1);
        end
        % recombine the spirals
        b = cat(1,bin, bout);
        %}

     
        %throw out ramp
        b = b(Nramp+1:end-Nramp,:,:);
        
        % 11.20.21:  throw out center points
        b(end/2-Ncenter/2:end/2+Ncenter/2, :,:) = [];
    
    

%         % save the center-crossings for diagnostics:
%         NCC = 20;
%         centers = b(end/2-NCC:end/2+NCC, :,:);
%         save center_xing.mat centers
        
        % rearrange the order of the echoes to match the acquisition order
        tmpk = [];
        for lf=1 : nleaves
            for sl=1 : nslices
                tmpk = [tmpk ; b(:,lf,sl)];
            end
        end
        
        tmpk(isnan(tmpk)) = 0;
        
        if SHOWPIX >2
            figure(1)
            subplot(221)
            plot((abs(tmpk(:)) / (max(abs(tmpk(:))))));
            title('concatenated echo train (reordered)')
            
            drawnow
        end
        
        
        % 3D gridding recon

        
        % interpolate to the cartesian grid:
        % adding weights to the data
        tmpk = tmpk(:).*weights(:);
        % 11.22.21: note the order change
        gkdata = griddata(ks(:,2), ks(:,1), ks(:,3), tmpk(:), kx, ky, kz);
       
        % alternative interpolation methods:
        % gkdata = interp3(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz, 'cubic');
        % F = scatteredInterpolant(ks(:,1), ks(:,2), ks(:,3), tmpk(:),'linear','none');
        % gkdata = F(kx,ky,kz)
        % dist = sqrt(kx.^2 + ky.^2 + kz.^2);
        % dist = 1;
        
        % where the interpolation returned NaN, we must put in zeros
        gkdata(isnan(gkdata)) = 0;
        gkdata(isinf(gkdata)) = 0;
            
             
        if ~doNUFFT
            % simple recon: do the FFT in 3D:
            % gkdata = gkdata ;
            
            im = fft3d(gkdata );
            
        else
            % 3D NUFFT recon

            im = nufft3d(-ks, tmpk(:), fov3, dim3, kdens);
            
        end
        
        % store in a buffer for parallel computing
        coilbuf(:,:,:,c) = (im);
        
        if SHOWPIX >1
            figure(1)
            subplot(223)
            lightbox(reshape(abs(im), dim3)); 
            title(sprintf('gridding recon - coil %d', c))
            colormap parula
            drawnow
            
            figure(1)
            subplot(224)
            lightbox(log(abs(gkdata)),[0 7],[]);
            title(sprintf('Log of re-grided k-data, coil %d', c))
            drawnow
        end
        
        if SHOWPIX > 2         
            figure(3)
            ov(ah,log(abs(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            %ov([],(angle(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            %set(gcf, 'Name', sprintf('Regridded k-space (log) - coil%d' , c));
            %set(gcf, 'Position', [600 0 0 0]);
            drawnow
        
            figure(2)
            ov(ah, abs(im), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            %ov([],(angle(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            set(gcf, 'Name', sprintf('Ortho views of reconned image, coil %d' , c));
            drawnow
        end

        
    end
    
    % now combine the coils (outside of the parallel loop)
    if ncoils==1
        im_all = coilbuf(:,:,:,1);
        
    else
        
        for c=1:ncoils
            
            % try forcing the mean phase of each coil images to agree.
            %         tmp = coilbuf(:,:,:,c);
            %         phscoil = angle(mean(tmp(:)));
            %         im_all = (im_all) + coilbuf(:,:,:,c) ./ exp(i * phscoil);
            
            im_all = abs(im_all) + i*abs(coilbuf(:,:,:,c));
            
        end
        
    end
    
    
    fprintf('\nframe %d ... Done! \n' , f);
    
    % flip image along z-axis:
    %im_all(:,:,[end:-1:1]) = im_all(:,:,[1:end]);
    im_all = im_all(end:-1:1);
    im_all = reshape(im_all, dim3);
    if SHOWPIX > 0
        %delete(figure(2)) ; delete(figure(3));
        figure(3);
        ov(ah,abs(im_all), round(dim3(1)/2 ), round(dim3(2)/2 ), round(dim3(3)/2 ), 0);
        set(gcf, 'Name', sprintf('Ortho views of reconned image %d', f));
        %fpos = get(gcf, 'Position');
        colormap bone
        drawnow
        
        %
        figure(4)
        gkdata = fft3d(im_all);
        ov(ah,log(abs(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
        set(gcf, 'Name', 'Ortho views of k space: log(FFT(image))' );
        drawnow
        %
        
    end
    imseries(:,f) =  im_all(:);
    
end



%%
h=define_nii_hdr();
h.dim = [4 dim3 nframes 0 0 0];
h.pixdim = [1 fov3./dim3 1 0 0 0];
h.datatype = 4;
h.bitsperpixel = 16;
% save myspace.mat
imseries = imseries/max(imseries(:))*((2^15)-1);
write_nii('timeseries_mag.nii', abs(imseries(:)), h,0);
write_nii('timeseries_angle.nii', angle(imseries(:))*1000, h,0);
if nframes>3
    h.dim(5) = 1;
    write_nii('timeseries_snr.nii', ...
        mean(imseries(:,3:end), 2)./ ...
        std(imseries(:,3:end),[],2), h,0);
end
return
%%
function dens = sampdensity(samplocs, N)
% function dens = sampdensity(samplocs, N)
%
% compute the sampling density of a set of points (samplocs)
% by summing the distance of the points N nearest neighbors
%
% samplocs are the coordinates of the samples usually N x 3
% N is the number of neighbors to evaluate
%
% the output is scaled so that the maximum density is 1
%

% first downsample the number of locations
F = 10;
N0 = [1:length(samplocs)]';
N1 = [1:F:length(samplocs)]';
samplocs = samplocs(1:F:end,:);

dens = zeros(size(samplocs,1), 1);

parfor n = 1:length(dens)
    distances = vecnorm((samplocs(n,:) - samplocs), 2, 2);
    distances = sort(distances);
    dens(n) = sum(distances(1:N));
end

dens = dens./max(dens(~isinf(dens)));
dens(isnan(dens)) = 0;
dens(isinf(dens)) = 1;

% upsample back to the original number
% we assume that the desity function is smooth
dens = abs(interp1(N1,dens,N0));
return


function dens = sampdensity_dist(samplocs, r)
% function dens = sampdensity(samplocs, r)
%
% compute the sampling density of a set of points (samplocs)
% by counting the number of points within a neighborhood of size r
%
% samplocs are the coordinates of the samples usually N x 3
% r is the radius of the neighborhood used to count samples
%
% the output is scaled so that the maximum density is 1
%

% first downsample the number of locations
F = 10;
N0 = [1:length(samplocs)]';
N1 = [1:F:length(samplocs)]';
samplocs = samplocs(1:F:end,:);

dens = zeros(size(samplocs,1), 1);

for n = 1:length(dens)
    for m =1:length(dens)
        if norm(samplocs(n, :) - samplocs(m,:)) < r
            dens(n) = dens(n) + 1;
        end
    end
end

dens = dens/max(dens);
% upsample back to the original number
% we assume that the desity function is smooth
dens = abs(interp1(N1,dens,N0));
return
%%
function out = phaseDetrend2(in, ph_ref)
% phase correction for stack of slpirals.  do not force phase at center of k space
% to be zero
% force phase at center of k-space to be ph_ref

m = abs(in);
p = angle(in);
p2 = p;
BUF = 20;

p = unwrap(p);

slope = (p(round(end/2+BUF)) - p(round(end/2-BUF))) / BUF/2;
lintrend = [0:length(p2)-1]*slope ;
p2 = p - lintrend';
p2 = p2 - p2(end/2) + ph_ref;

%{
plot(p, 'b')
hold on
plot(p2, 'r')
hold off
drawnow
%}
out = m.* exp(-1i*p2);
return

%%
function out = phaseDetrend(in)
% function out = phaseDetrend(in)
% takes a complex signal and subtracts a linear trend from the phase component
% we then force the center of the line to go through zero

m = abs(in);
p = angle(in);
p2 = p;
BUF = 20;
%{
% (correction 1) use the center points to establish a line slope
%  preserve the phase in the mid point
p = unwrap(p);
slope = (p(end/2+BUF)-p(end/2-BUF)) / BUF/2;
lintrend = [0:length(p)-1]*slope;
lintrend = lintrend - lintrend(end/2);  % zero phase in the mid point.
p2 = p - lintrend' ;
%}
%{
% (correction 2)linear detrend the whole echo
p = unwrap(p);
p2 = detrend(p,1) ;
p2 = p2-p2(end/2);
%}
% (correction 3) use the end points to establish a line slope
%{
slope = (p(end)-p(1)) / length(p);
lintrend = p(end/2) + [0:length(p)-1]*slope;
p2 = p -lintrend' ;
%}

%
% (correction 4 *) phase shift the whole echo so that the center is always 0
p = unwrap(p);

phase_cor = mean(p2(round(end/2-5 : end/2+5)));
p2 = p - phase_cor;

out = m.* exp(1i*p2);
%{
plot(p, 'b')
hold on
plot(p2, 'r')
hold off
pause
%}

%{
% just average together all the center points
tmp = mean(in(end/2-10:end/2+10));
out = in;
out(end/2-10:end/2+10) =tmp;
%}


%{
% (correction 5) Real and imaginary: use the center points to establish a line slope
%  preserve the value in the mid point
re = real(in);
im = imag(in);

slope = (re(BUF)-re(end-BUF)) / (length(in) - 2*BUF);
lintrend = [0:length(re)-1]*slope;
lintrend = lintrend - lintrend(round(end/2));  % zero phase in the mid point.
re2 = re - lintrend' ;

slope = (im(BUF)-im(end-BUF)) / (length(in) - 2*BUF);
lintrend = [0:length(re)-1]*slope;
lintrend = lintrend - lintrend(round(end/2));  % zero phase in the mid point.
im2 = im - lintrend' ;
%
plot(re, 'b')
hold on
plot(re2, 'r')
hold off
grid on
pause
%

out = complex(re2, im2);
%}
%save phase_trend.txt slope -ascii -append

return


%{
function area = voronoidens (kxyz)
% function area = voronoidens(kxyz);
% input: kxyz is k-space trajectories
% output: area of cells for each point

% create complex k-space locations and find unique points

[C, ia, ic] = unique(kxyz);
nunique_pts = numel(C);

% returns vertices and cells of
% voronoi diagram
[V,C] = voronoin(kxyz);

% unpack cell array, compute area of each ploygon
areas_unique = zeros(nunique_pts,1);
for jj = 1:length(C)
    x = V(C{jj},1);
    y = V(C{jj},2);
    z = V(C{jj},3);
    
    x(isinf(x)) = 1e5;
    y(isinf(y)) = 1e5;
    z(isinf(z)) = 1e5;
    
    [k, v] = convhull(x,y,z);
    areas_unique(jj) = v;
end

area = areas_unique(ic);

return
%}

%%  extra code for debugging
        
        %%  2D test code
%NUFFT  recon:
%{
    bufsize = size(tmpk,1);
    for s=1:nslices
        kxy = ks(bufsize*(s-1)+1:bufsize*s, 1:2);
        kdat = tmpk(:,s);
        xcp(:,:,s) = nufft2d(kxy, kdat, fov3(1:2), dim3(1:2));
    end
    subplot(224)
    imagesc(abs(xcp(:,:,1)'))
    title('2D NUFFT recon')
%}

%  2D gridding recon:
%{
    [kx , ky, kz ]= meshgrid( ...
        linspace(-Kxymax,Kxymax,MYDIM), ...
        linspace(-Kxymax,Kxymax,MYDIM), ...
        linspace(-Kzmax,Kzmax,nslices));
    
    for s=1:nslices
        F = scatteredInterpolant(ks(:,1), ks(:,2), tmpk(:,s),'linear','none');
        gkdata = F(kx(:,:,s),ky(:,:,s));
        gkdata(isnan(gkdata)) = 0;
        xcp(:,:,s) = fftshift(fft2(fftshift(gkdata)));
    end
    subplot(224)
    lightbox(log(abs(gkdata)),[],[]); title('Log k-space image (re-grided)')
    
    subplot(222)
    lightbox(fliplr(abs(xcp)))
    title('2D gridding recon')
    return
%}

%% code for trajectory mapping:
%{
        if doTraj
            figure(1)
            
            buf = raw(f,:);
            buf = reshape(buf,ndat,nslices);
            buf = buf(:,nslices/2);
            %buf = buf(:,3);
            
            if n==1 % this is a regular frame
                whole = buf;
                
            end
            
            if n==2 % this is the reference frame: no grads.
                ref = buf;
                subplot(222)
                hold off
                plot(phase(whole)/gamma_rad/dt2/1e3)
                hold on
                plot(unwrap(angle(ref)))
                hold off
                title('frames 1 and 2: signal phases')
                legend('both grads on', 'no grads (ref)')
            end
            
            if n==3 % x gradient only
                buf = buf ./ref;
                kx_measured = unwrap(angle(buf))/delta_x/2/pi;
                subplot(223)
                hold off
                plot(kx_measured)
                hold on
                plot(ks(1:ndat,1));
                hold off
                legend('kx measured',  'kx nominal')
                title(sprintf('Kx (frame 3)'))
                
            end
            if n==4
                buf = buf ./ref;
                ky_measured = unwrap(angle(buf))/delta_y/2/pi;
                subplot(224)
                hold off
                plot(ky_measured)
                hold on
                plot(ks(1:ndat, 2));
                legend(  'ky measured', 'ky nominal')
                title(sprintf('Ky (frame 4)'))
                hold off
                
            end
        end
%}
