function [tmpk, im_all] = spiral3drecon_201119(myPfile,doRots, Ndel,  doNUFFT)
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
scaleFactor = 1;
doT2filter = 0;
T2 = 0.150;
Nramp = 0;         % Number of ramp points to throw out
SHOWPIX = 3;
doPhaseCyle = 1;   % CPMG phase cycling
doPhaseDetrend = 0;% linear phase correction for each echo
useKrot = 1;       %   useKrot:  read in krots.txt and use it to calculate the rotations from interleaves v. slice encodes
doGoldenAngle = 1; %   doGoldenAngle: rotations done by Golden Angle (137.5) vs. 180/n.slices?
!rm phis_file.mat

if nargin==0
    pfiles = dir('P*.7');
    myPfile = pfiles(1).name;
    doNUFFT = 1;
    doRots = 1;
    doGoldenAngle = 1;
    Ndel = -3;
    useKrot = 0;
    
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

%%%% for debugging:


if doRots
    fov3 = [fov fov fov];
    dim3 = [MYDIM MYDIM MYDIM];
end
dim3 = dim3+1;  % the center of k-space is important, so make sure the dimensions allow it

% make some headers for output
h=define_nii_hdr();
h.dim = [4 dim3 nframes 0 0 0];
h.pixdim = [1 fov3./dim3 1 0 0 0];
h.datatype = 4;
h.bitsperpixel = 16;
ah = nii2avw_hdr(h);

if SHOWPIX >2
    delete(figure(1))
    figure(1)
    subplot(221)
    plot(squeeze((abs(raw(1,:,:)))));
    title('Echo train: Frame 1, all coils');
    fpos = get(gcf,'Position');
    set(gcf,'Position', fpos + [0 600 0 0]);
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
Ndel2= 0;


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

if doRots  % SERIOS readout: rotating spiral pancakes about x axis
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
        krot = krot(:, 5:end);
        krot = krot/32767;
        
    else
        % calculate rotation angles on the fly, 
        % assuming GoldenAngle rotations about x-axis and (slices)
        % evenly spaced rotations about z-axis (leaves)
        angle_list = zeros(nslices*nleaves,2);
        zrot_angle = pi/nleaves;
        count = 1;
        for n=0:nleaves-1
            for m=0:nslices-1
                angle_list(count,1) = deg2rad(137.5)*m;
                angle_list(count,2) = zrot_angle*n;
                count = count +1;
            end
        end
    end
    
    ks_all = [];
    echolist = [];
    idx=0;
    
    
    for s = 0:nslices-1
        
        ks_plate = ks';
        
        for leafn = 0:nleaves-1
            
            glen = size(ks_plate,1);
            idx = leafn*nslices + s +1 ;
            if 0
                % use the rot matrices from file ... this doesn't work
                % (Rx * Rz is incorrect, but it's what the scanner wants)
                % is scaling the problem??
                row = krot(idx,:);
                rot = reshape(row, 3,3) ;
                tmpks = rot*ks_plate';
                
            else
                % use the angles from the file
                % and generate rot matrices from those
                % rotation angles are in the first two columns of the file
                xi = angle_list(idx,1);
                phi = angle_list(idx,2);
                
                % 8.6.2020: the rotation matrix is doing (-Rx)*Rz
                % xi needs to be negative in order for rotations to agree
                % with angles
                sxi = sin(xi); cxi = cos(xi);
                sphi = sin(phi); cphi = cos(phi);
                
                Rz = [cphi -sphi 0;
                    sphi cphi 0;
                    0 0 1];
                
                Rx = [1  0   0 ;
                    0 cxi sxi ;
                    0 -sxi cxi];
                
                tmpks = Rz * ks_plate';
                tmpks = Rx * tmpks;
                
                echolist = [echolist; s, leafn , xi, phi];
            end
            
            ks_all = [ks_all; tmpks'];
            
            if SHOWPIX > 3
                figure(1)
                subplot(2,2,2)
                plot3(tmpks(1,:), tmpks(2,:), tmpks(3,:));
                title(sprintf('Order : slice= %d, leaf= %d', s, leafn));
                axis([-1 1 -1 1 -1 1]*Kxymax*3)
                drawnow
                hold off
            end
        end
    end
    ks = ks_all;
    
    save phis_file.mat angle_list echolist
    
    
    
    
else  % stack of spirals
    
    % calculate the kz components for each echo
    Kzmax = 1/(slthick)/2;
   
    % the kz steps order and corresponding fraction of Gz max are
    % specified in a file
    kzsteps = load('kzsteps.txt');
    kzfraction = kzsteps(:,3);      % faction of kz space for each platter
    kzorder = kzsteps(:,2) +1;      % order of the platters during readout  (center -> out)
%     for n=1:numel(kzfraction)
%         kzfraction(kzorder(n)) = kzsteps(n,3);      % faction of kz space for each platter
%     end
 
 
     kzfraction = -kzfraction;
     
    % calculate kx, ky kz locations for each platter:
    % kx and ky are always the same, but kz needs to 
    % be scaled for every platter
    
    % copy the xy trajectory to all platters
    tmp = repmat(ks,  1, nslices*nleaves);
    kxs = tmp(1,:);
    kys = tmp(2,:);
    kzs = tmp(3,:);
    
    
    % now re-scale:
    fin = 0;
    for n=1:nslices
        
        % if  more than one interleave per platter,
        % we need to rotate about the z-axis for each
        % interleaf.
        for m=1:nleaves
            phi = (m-1)*pi/nleaves;
            sphi = sin(phi); cphi = cos(phi);
            
            Rz = [cphi -sphi 0;
                sphi cphi 0;
                0 0 1];
            Rtmp = Rz * tmp;            
            
            beg = fin+1;
            fin = beg + (ndat-2*Nramp)-1;
            
            % [beg fin (fin-beg) ndat length(ks) length(tmp) length(kxs)]
            
            kxs(beg:fin) = Rtmp(1,beg:fin);
            kys(beg:fin) = Rtmp(2,beg:fin);
            
            kzs(beg:fin) = kzs(beg:fin) * kzfraction(n);

            if SHOWPIX >3
                figure(1)
                subplot(2,2,2)
                plot3(kxs(beg:fin), kys(beg:fin), kzs(beg:fin));
                axis(Kxymax*[-1 1 -1 1 -1 1]*Kxymax);
                hold off
                pause(0.1)
                drawnow
            end
        end
    end
    
    % rejoin the coordinates of the trajectory into a single matrix
    ks = [kxs' kys' kzs'];
    
    % 2D test
    % kzs(beg:fin) = 0;
    
end

if SHOWPIX >2
    figure(1)
    subplot(2,2,2)
    plot3(ks(:,1), ks(:,2), ks(:,3));
    axis([-1 1 -1 1 -1 1]*Kxymax*3)
    hold off
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

MaxFrame = nframes
if SHOWPIX==3
    MaxFrame = 1;
end

raw = raw * scaleFactor;

for f=1 :MaxFrame
    
    im_all = 0;
    coilbuf = zeros(dim3(1), dim3(2), dim3(3), ncoils);
    
    parfor c = 1:ncoils
        fprintf('\nReconning frame %d , coil %d, delay %d...\n', f, c, Ndel);
        
        % reshape, filter and adjust the data:
        
        tmpk = raw(f,:,c);   % new dims:  Nslices*ndat*nleaves x 1
        tmpk = reshape(tmpk, ndat*nleaves, nslices);
        
      
        %throw out ramp
        b = reshape(tmpk,ndat,nleaves,nslices);
        b = circshift(b, [Ndel 0 0]);
        b = b(Nramp+1:end-Nramp,:,:);
        
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
        end
        
        % Remove linear phase trend in the data
        if doPhaseDetrend
            for ss=1:nslices
                b(:,:,ss) = phaseDetrend(b(:,:,ss));
            end
            
        end
        
        % add a T2 filter
        if doT2filter
            if ~doRots
                % the reference is the first few points ot the first echo
                % rescale other echoes so that the center of k-space has
                % the same amplitude as the first echo
                kcenter0 = mean(abs(b(1:10, 1,1)));
                for ss = 1:nslices
                    for ll=1:nleaves
                        kcenter = mean(abs(b(1:10, ll,ss)));
                        b(:,ll,ss) = b(:,ll,ss) * kcenter0 / kcenter;
                    end
                end
                
            else
                % the reference is a few points at the center ot the first echo
                % rescale other echoes so that the center of k-space has
                % the same amplitude as the first echo
                kcenter0 = mean(abs(b(end/2-5:end/2+5, 1,1)));  
                for ss = 1:nslices
                    for ll=1:nleaves
                        kcenter = mean(abs(b(end/2-5:end/2+5, ll,ss)));
                        b(:,ll,ss) = b(:,ll,ss) * kcenter0 / kcenter;
                    end
                end
            end
        end
        
        tmpk = reshape(b,(ndat-2*Nramp)*nleaves,nslices);
        
        tmpk(isnan(tmpk)) = 0;
        
        if SHOWPIX >2
            figure(1)
            subplot(221)
            plot(ks(:,2)/Kxymax)
            plot(ks(:,3)/Kzmax)
            hold on
            plot((abs(tmpk(:)) / (max(abs(tmpk(:))))));
            hold off
            title('echo train (overlaid on kx wave)')
            drawnow
        end
        
        
        % 3D gridding recon

        
        % interpolate to the cartesian grid:
        gkdata = griddata(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz);
        
        % alternative interpolation methods:
        % gkdata = interp3(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz, 'cubic');
        % F = scatteredInterpolant(ks(:,1), ks(:,2), ks(:,3), tmpk(:),'linear','none');
        % gkdata = F(kx,ky,kz)
        % dist = sqrt(kx.^2 + ky.^2 + kz.^2);
        % dist = 1;
        
        % where the interpolation returned NaN, we must put in zeros
        gkdata(isnan(gkdata)) = 0;
        gkdata(isinf(gkdata)) = 0;
        
        if SHOWPIX >2
            subplot(224)
            lightbox(log(abs(gkdata)),[],[]);
            title(sprintf('Log of re-grided k-data, coil %d', c))
            
        end
        
        if SHOWPIX > 1
            
            figure(3)
            ov(ah,log(abs(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            %ov([],(angle(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
            set(gcf, 'Name', sprintf('Ortho views of regridded k-space (log) %d' , f));
            set(gcf, 'Position', fpos + [600 0 0 0]);
            drawnow
        end
        
        if ~doNUFFT
            % simple recon: do the FFT in 3D:
            % adding weights to the data
            % gkdata = gkdata; % .* (1+dist).^0.5;
            
            im = fft3d(gkdata );
            
        else
            % 3D NUFFT recon
            % make weighting function proportional to distance to center of
            % kspace
            if doRots
                kdens= sqrt(ks(:,1).^2 + ks(:,2).^2 + ks(:,3).^2 );
                kdens = kdens/(max(kdens));
                kdensfac =  1.5;
                
                %
                % let's explore looking at the trajectory speed for density
                % compensation
                dks = diff(ks,2);
                kdens2 = sqrt(dks(:,1).^2 + dks(:,2).^2 + dks(:,3).^2 );
                kdens2 = [kdens2(1,:) ; kdens2; kdens2(end,:)];
                kdens2 = kdens2/(max(kdens2));
                kdens = kdens2 + kdens;
                kdensfac = 3;
                %}
            else
                kdens= sqrt(ks(:,1).^2 + ks(:,2).^2 );
                kdensfac = 0.9;
            end
            kdens = kdens .^ kdensfac;
            kdens = 1e-6 + kdens/max(kdens);
            
            im = nufft3d(ks, tmpk(:), fov3, dim3, kdens);
            
        end
        
        % store in a buffer for parallel computing
        coilbuf(:,:,:,c) = im;
    end
    
    % now combine the coils (outside of the parallel loop)
    for c=1:ncoils
        %im_all = abs(im_all) + i*abs(coilbuf(:,:,:,c));
        
        % try forcing the mean phase of each coil images to agree.
        tmp = coilbuf(:,:,:,c);
        phscoil = angle(mean(tmp(:)));
        im_all = (im_all) + coilbuf(:,:,:,c) ./ exp(i * phscoil);
        
        if SHOWPIX >2
            subplot(223)
            lightbox(abs(im_all)); title('gridding recon')
            colormap parula
            drawnow
        end
        
    end

    fprintf('\nframe %d ... Done! \n' , f);
    
    % flip image along z-axis:
    im_all(:,:,[end:-1:1]) = im_all(:,:,[1:end]);
    
    if SHOWPIX > 0
            
        delete(figure(2)) ; delete(figure(3));
        figure(2);
        ov(ah,abs(im_all), round(dim3(1)/2 ), round(dim3(2)/2 ), round(dim3(3)/2 ), 0);
        set(gcf, 'Name', sprintf('Ortho views of reconned image %d', f));
        fpos = get(gcf, 'Position');
        colormap bone
        drawnow
    end

    imseries(:,f) =  im_all(:);
    
end
%%
h=define_nii_hdr();
h.dim = [4 dim3 nframes 0 0 0];
h.pixdim = [1 fov3./dim3 1 0 0 0];
h.datatype = 4;
h.bitsperpixel = 16;
save myspace.mat
write_nii('timeseries_mag.nii', abs(imseries(:))*1000, h,0);
write_nii('timeseries_angle.nii', angle(imseries(:))*1000, h,0);

return
%%
function out = phaseDetrend(in)
% function out = phaseDetrend(in)
% takes a complex signal and subtracts a linear trend from the phase component
% the linear trend is calculated from the first and last points in the data
% assuming that they are supposed to bt the same

m = abs(in);
p = angle(in);

slope = (p(end)-p(1)) / length(p);
lintrend = p(1) + [0:length(p)-1]*slope;

p2 = p -lintrend' ;

% p2 = detrend(p,3) ;

out = m .* exp(1i*p2);
save phase_trend.txt slope -ascii -append

return

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
