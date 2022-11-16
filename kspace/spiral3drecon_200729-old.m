function [echo, im_all] = spiral3drecon_200729(doRots, useKrot, Ndel, doGoldenAngle)
% function im_all = spiral3drecon_200729(doRots, useKrot, Ndel, doGoldenAngle)
%
% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot - USES rotation marices sotred in text files.
% 3 - reconstructs the image by Gridding and by Jeff's NUFFT
%

dt = 4e-6; %seconds
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

doT2filer = 0;
doNUFFT = 0;
Nramp = 0;
SHOWPIX = 1;

T2 = 0.150;

pfiles = dir('P*.7');
myPfile = pfiles(1).name;

if nargin==0
    doRots = 1;
    doGoldenAngle = 0;
    doTraj= 0;
    Ndel = 0;
    useKrot = 1;
    
    delta_x = 7; % plane offset for Trajectory measurement
    delta_y = 7; % plane offset for Trajectory measurement
end
fprintf('\nBegin recon of 3D spiral trajectory for %s, delay=%d ...', myPfile, Ndel)
%run('~hernan/matlab/irt/setup')

% load the raw data:
% Nframes x Nslices*ndat x Ncoils
[raw scaninfo ] = read_raw_3d(myPfile ,0 );


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

if doRots
    fov3 = [fov fov fov];
    dim3 = [MYDIM MYDIM MYDIM];
end
dim3 = dim3+1;  % the center of k-space is important, so make sure the dimensions allow it

if SHOWPIX
    figure(1)
    subplot(221)
    plot(squeeze((abs(raw(1,:,:)))));
    title('Echo train: Frame 1, all coils');
    drawnow
end

echo = squeeze((abs(raw(1,:,:))));



% Now reshape into a 3D matrix
% kdata = reshape(kdata,  ndat, nslices, nframes);

fprintf('\nLoading grads...');
g = load('grad.txt');
if useKrot
    % g = [g(:,2) , g(:,1) , g(:,3)];
    % the scanner frame of reference has the x and y channels flipped
    % verified in 2D image that this in needed
end
glen = size(g,1);
aqwindow = dt*glen;
dt2 = aqwindow/ndat;
%%


cnt = 1;

Ndel2= 0;


% use the k-space trajectory calculated from the gradients,
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;  % this is in cm-1 units

% tranpose to make things work
k2 =  k2';

Kxymax = max(abs(k2'));
Kxymax = max(Kxymax(1:2));
Kxymax = MYDIM/fov/2;  % Kmax
Kzmax = Kxymax;

% must resample the trajectory to match the data acquisition rate
ks = zeros(3,ndat);
fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
for n=1:3
    tmp = k2(n,:);
    %ks(n,:) = tmp(1:ndat);
    ks(n,:) = interp1(linspace(0,1,length(tmp)) , tmp, linspace(0,1,ndat));
end
% use interp

% throw out ramp
ks = ks(:, Nramp+1:end-Nramp) ;

if doRots  % rotating spiral pancakes about x axis
    
    
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
        % the last 9 columns are the actual rotation matrix
        
        krot = load('krotations.txt');
        angle_list = krot(:,3:4);
        krot = krot(:, 5:end);
        krot = krot/32767;
        
    else
        % calculate rotation angles on the fly:
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
            %             figure(1)
            %             subplot(2,2,2)
            %             plot3(tmpks(1,:), tmpks(2,:), tmpks(3,:));
            %             title(sprintf('Order : slice= %d, leaf= %d', s, leafn));
            %             axis([-1 1 -1 1 -1 1]*Kxymax)
            %             drawnow
            
            hold on
        end
    end
    hold off
    ks = ks_all;
    
    save phis_file.mat angle_list echolist
    if SHOWPIX
        figure(1)
        subplot(2,2,2)
        plot3(ks(:,1), ks(:,2), ks(:,3));
        title(sprintf('All echoes: slices= %d, leaves= %d', s, leafn));
        axis([-1 1 -1 1 -1 1]*Kxymax)
        hold off
    end
    
    
    
else  % stack of spirals
    
    % the kz steps order and corresponding fraction os Gz max are
    % specified in a file
    kzsteps = load('kzsteps.txt');
    kzfraction = kzsteps(:,3);
    kzorder = kzsteps(:,2);
    ztmp = zeros(size(kzfraction));
    for n=1:nslices
        p = kzorder(n)+1;
        ztmp(p) = kzfraction(n);
    end
    kzfraction = ztmp';
    
    % calculate kx, ky kz locations for
    % kx and ky are always the same, but kz increases monotonically
    % this was not recorded in the grad file
    % copy the xy trajectory to all platters
    tmp = repmat(ks,  1, nslices*nleaves);
    kxs = tmp(1,:);
    kys = tmp(2,:);
    kzs = tmp(3,:);
    
    % calculate the kz components for each echo
    Kzmax = 1/2/(slthick);
    
    for m=1:nleaves
        phi = m*pi/nleaves;
        sphi = sin(phi); cphi = cos(phi);
        
        Rz = [cphi -sphi 0;
            sphi cphi 0;
            0 0 1];
        Rtmp = Rz * tmp;
        
        for n=1:nslices
            beg = (m-1)*(ndat-2*Nramp)*(n-1)+1;
            fin = (m-1)*(ndat-2*Nramp)*n;
            kxs(beg:fin) = Rtmp(1,beg:fin);
            kys(beg:fin)= Rtmp(2,beg:fin);
            kzs(beg:fin) = kzs(beg:fin)*kzfraction(n);
            
            if SHOWPIX
                figure(1)
                subplot(2,2,2)
                plot3(kxs(beg:fin), kys(beg:fin), kzs(beg:fin));
                axis(Kxymax*[-1 1 -1 1 -1 1]);
                hold on
                drawnow
            end
        end
    end
    
    
    % 2D test
    % kzs(beg:fin) = 0;
end
hold off
    % rejoin the components into a single matrix
    ks = [kxs' kys' kzs'];
    
    


% the dimensions of raw are: Nframes x Nslices*ndat x Ncoils
imseries = zeros(dim3(1)*dim3(2)*dim3(3),nframes);

for f=1:1
    im_all = 0;
    for c = 1:ncoils
        fprintf('\nReconning frame %d , coil %d, delay %d...\n', f, c, Ndel);
        
        tmpk = raw(f,:,c);   % new dims:  Nslices*ndat*nleaves x 1
        tmpk = reshape(tmpk, ndat*nleaves, nslices);
        
        %throw out ramp
        b = reshape(tmpk,ndat,nleaves,nslices);
        b = circshift(b, [Ndel 0 0]);
        b = b(Nramp+1:end-Nramp,:,:);
        
        %
        % add a T2 filter
        if doT2filer
            if ~doRots
                for ss = 1:nslices
                    b(:,:,ss) = b(:,:,ss)*exp(2*abs(ss-nslices/2-1)*15/125);
                end
            else
                for ss = 1:nslices
                    b(:,:,ss) = b(:,:,ss)*exp((ss-1)*15/125);
                end
            end
        end
        %}
        tmpk = reshape(b,(ndat-2*Nramp)*nleaves,nslices);
        
        tmpk(isnan(tmpk)) = 0;
        
        if SHOWPIX
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
        
        %
        
        % 3D gridding recon
        [kx , ky, kz ]= meshgrid( ...
            linspace(-Kxymax,Kxymax,dim3(1)), ...
            linspace(-Kxymax,Kxymax,dim3(2)), ...
            linspace(-Kzmax,Kzmax,dim3(3)));
        
        dist = sqrt(kx.^2 + ky.^2 + kz.^2);
        %dist = 1;
        %F = scatteredInterpolant(ks(:,1), ks(:,2), ks(:,3), tmpk(:),'linear','none');
        %gkdata = F(kx,ky,kz);
        
        % alternative interpolation
        gkdata = griddata(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz);
        % gkdata = interp3(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz, 'cubic');
        
        % where the interpolation returned NaN, we must put in zeros
        gkdata(isnan(gkdata)) = 0;
        gkdata(isinf(gkdata)) = 0;
        
        if SHOWPIX
            subplot(224)
            lightbox(log(abs(gkdata)),[0 7],[]);
            title(sprintf('Log of re-grided k-data, coil %d', c))
        end
        
        if ~doNUFFT
            %  do the FFT :
            gkdata = gkdata; % .* (1+dist).^0.5;
            im = fft3d((gkdata ));
        else
            % 3D NUFFT recon
            im = nufft3d(ks, (tmpk(:)), fov3, dim3);
        end
        
        % combine the coils:
        im_all = abs(im_all) + i*abs(im);
        
        if SHOWPIX
            subplot(223)
            lightbox(abs(im_all)); title('gridding recon')
            colormap parula
            drawnow
        end
        
    end
    %%
    
    %%
    fprintf('\nframe %d ... Done! \n' , f);
    
    % flip image along z-axis:
    im_all(:,:,[end:-1:1]) = im_all(:,:,[1:end])
    
    figure(2);
    ov([],abs(im_all), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
    title(sprintf('Ortho views of reconned image %d', f))
    
    figure(3)
    ov([],log(abs(gkdata)), round(dim3(1)/2), round(dim3(2)/2), round(dim3(3)/2), 0);
    title(sprintf('Ortho views of regridded k-space (log) %d' , f))
    
    imseries(:,f) =  im_all(:);
    
    
end
%%
h=define_nii_hdr();
h.dim = [4 dim3 nframes 0 0 0];
h.pixdim = [1 fov3./dim3 1 0 0 0];
h.datatype = 4;
h.bitsperpixel = 16;
save myspace.mat
write_nii('timeseries_mag.nii', abs(imseries(:))/100, h,0);
write_nii('timeseries_angle.nii', angle(imseries(:))*1000, h,0);

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
