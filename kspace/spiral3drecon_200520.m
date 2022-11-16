function im_all = spiral3drecon_200520(doRots, Ndel, doGoldenAngle, doTraj)
% function im_all = spiral3drecon_200520(doRots, Ndel, doGoldenAngle, doTraj)
%
% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - reconstructs the image by Gridding and by Jeff's NUFFT
%

dt = 4e-6; %seconds
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

MYDIM = 64;


pfiles = dir('P*.7');
myPfile = pfiles(1).name;

if nargin==0
    doRots = 1;
    doGoldenAngle = 0;
    doTraj= 0;
    Ndel = -1;
    
    delta_x = 7; % plane offset for Trajectory measurement
    delta_y = 7; % plane offset for Trajectory measurement
end

fprintf('\nBegin CP NUFFT recon of spiral trajectory for %s, delay=%d', myPfile, Ndel)
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
    %fov3 = [fov fov fov];
    dim3 = [MYDIM MYDIM MYDIM];
end

figure(1)
subplot(221)
plot(squeeze((abs(raw(1,:,:)))));
title('Echo train: Frame 1, all coils');
drawnow



% Now reshape into a 3D matrix
% kdata = reshape(kdata,  ndat, nslices, nframes);

fprintf('\nLoading grads...');
g = load('grad.txt');
glen = size(g,1);
aqwindow = dt*glen;
dt2 = aqwindow/ndat;
%%
Nramp = 10;
slord = [];
m=1;
for n=1:nslices/2
    slord(m) =  nslices/2 + (n);
    slord(m+1) =  nslices/2 +1-(n);
    m=m+2;
end

imseries = [];
cnt = 1;

Ndel2= 0;
figure(1)


% use the k-space trajectory calculated from the gradients,
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;  % this is in cm-1 units

% tranpose to make things work
k2 =  k2';

Kxymax = max(abs(k2'));
Kxymax = max(Kxymax(1:2));
Kxymax = MYDIM/fov/2;  % Kmax

% must resample the trajectory to match the data acquisition rate
ks = zeros(3,ndat);
fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
for n=1:3
    tmp = k2(n,:);
    % ks(n,:) = interp1(linspace(0,1,glen), tmp, linspace(0,1,ndat));
    ks(n,:) = tmp(1:ndat);
end

% throw out ramp
ks = ks(:, Nramp+1:end-Nramp) ;


% if multiple interleaves,
% generate the rotated version and then concatenate the rotated trajectories
if nleaves >1
    rotang = pi/nleaves;
    ks_all = [];
    for n=0:nleaves-1
        c = cos(rotang*n); s=sin(rotang*n);
        R = [c s 0; 
            -s c 0; 
            0 0 1];
        kstmp = R*ks;
        
        ks_all = [ks_all kstmp];
    end
    ks = ks_all;  % (ks now contains all the interleaves for a single slice) 
end

% calculate kx, ky kz locations for all the stack of spirals
% kx and ky are always the same, but kz increases monotonically
% this was not recorded in the grad file

subplot(222)
plot(ks(1,:), ks(2,:));
title('Trajectory k-space traj. from grad file ');

if doRots==0
    % copy the xy trajectory to all platters
    tmp = repmat(ks,  1, nslices);
    kxs = tmp(1,:);
    kys = tmp(2,:);
    kzs = tmp(3,:);
    
    % calculate the kz components
    Kzmax = 1/(2*slthick);
    kzz = 1*linspace(-Kzmax, Kzmax,nslices);
    for n=1:nslices
        beg = nleaves*(ndat-2*Nramp)*(n-1)+1;
        fin = nleaves*(ndat-2*Nramp)*n;
        kzs(beg:fin) = kzz(n);
        % 2D test
        % kzs(beg:fin) = 0;
    end
    % rejoin the components into a single matrix
    ks = [kxs' kys' kzs'];
else
    tmpks = [];
    slorder = [0:nslices-1];
    %slorder = 1:nslices;
    xrot_angles = pi*slorder/nslices;
    
    if (doGoldenAngle)
        xrot_angles = slorder*137.5*(pi/180);
    end
    
    % rotating spiral: isotropic k space sampling
    Kzmax = Kxymax;
    
    for xi= xrot_angles
        cxi = cos(xi); sxi = sin(xi);
        
        Rx =[ 1  0   0 ;
            0 cxi  -sxi ;
            0 sxi cxi];
        
        tmp =  Rx*ks;
        plot3(tmp(1,:), tmp(2,:), tmp(3,:));
        
        axis(Kxymax*[-1 1 -1 1 -1 1]);     
        %hold on
        drawnow
        tmpks = [tmpks tmp];
        
    end
    hold off
    ks = tmpks';
    
end
% ks is ordered as :
%     slice 1
%         leaf 1
%         leaf 2 ...
%     slice 2
%         leaf 1
%         leaf 2 ...
            

%plot3(ks(:,1), ks(:,2), ks(:,3));


% including a filter to un-do the T2 effect
T2 = 150;


%}

% raw is Nframes x Nslices*ndat x Ncoils
for f=1 %:nframes
    im_all = 0;
    for c = 1:ncoils
        fprintf('\nReconning frame %d , coil %d, delay %d...\n', f, c, Ndel);
        
        tmpk = raw(f,:,c);   % new dims:  Nslices*ndat*nleaves x 1
        tmpk = reshape(tmpk,ndat*nleaves,nslices);
        %throw out ramp
        %
        b = reshape(tmpk,ndat,nleaves,nslices);
        b = circshift(b, [Ndel 0 0]);

        b = b(Nramp+1:end-Nramp,:,:);
        %b(:,3,:)= 0;
        %{
        % add a T2 filter
        for ss=1:nslices
            b(:,:,ss) = b(:,:,ss)*(exp(ss*30/T2));
        end
        %}
        tmpk = reshape(b,(ndat-2*Nramp)*nleaves,nslices);
        
        tmpk(isnan(tmpk)) = 0;
        
        figure(1)
        
        subplot(221)
        plot(ks(:,2)/Kxymax)
        hold on
        plot((abs(tmpk(:))) / (max((abs(tmpk(:))))));
        hold off
        title('echo train (overlaid on kx wave)')
        drawnow
        
        
        %%
        % 2D NUFFT  recon:
        %                 bufsize = size(tmpk,1);
        %                 for s=1:nslices
        %                     kxy = ks(bufsize*(s-1)+1:bufsize*s, 1:2);
        %                     kdat = tmpk(:,s);
        %                     xcp(:,:,s) = nufft2d(kxy, kdat, fov3(1:2), dim3(1:2));
        %                 end
        %                 subplot(224)
        %                 imagesc(abs(xcp(:,:,1)'))
        %                 title('2D NUFFT recon')
        %
        %
        %                 %  2D gridding recon:
        %                 [kx , ky, kz ]= meshgrid( ...
        %                     linspace(-Kxymax,Kxymax,MYDIM), ...
        %                     linspace(-Kxymax,Kxymax,MYDIM), ...
        %                     linspace(-Kzmax,Kzmax,nslices));
        %
        %                 for s=1:nslices
        %                     F = scatteredInterpolant(ks(1:end/2,1), ks(1:end/2,2), tmpk(:,s),'linear','none');
        %                     gkdata = F(kx(:,:,s),ky(:,:,s));
        %                     gkdata(isnan(gkdata)) = 0;
        %                     xcp(:,:,s) = fftshift(fft2(fftshift(gkdata)));
        %                 end
        %                 subplot(221)
        %                 lightbox(log(abs(gkdata)),[],[]); title('Log k-space image (re-grided)')
        %                 %
        %                 subplot(223)
        %                 imagesc(fliplr(abs(xcp(:,:,2))))
        %                 title('2D gridding recon')
        
        % 3D NUFFT recon
        %{
        xcp = nufft3d(ks, (tmpk(:)), fov3, dim3);
        subplot(224)
        lightbox(abs(xcp));
        title('3D NUFFT recon')
        %}
          
        %Kzmax = 1/(slthick);

        % 3D gridding recon
        [kx , ky, kz ]= meshgrid( ...
            linspace(-Kxymax,Kxymax,dim3(1)+1), ...
            linspace(-Kxymax,Kxymax,dim3(2)+1), ...
            linspace(-Kzmax,Kzmax,dim3(3)+1));
        
        dist = sqrt(kx.^2 + ky.^2 + kz.^2);
        
%         F = scatteredInterpolant(ks(:,1), ks(:,2), ks(:,3), tmpk(:),'linear','none');
%         gkdata = F(kx,ky,kz);
        
        % alternative interpolation
        gkdata = griddata(ks(:,1), ks(:,2), ks(:,3), tmpk(:), kx, ky, kz);
        
        % where the interpolation returned NaN, we must put in zeros
        gkdata(isnan(gkdata)) = 0;
        gkdata(isinf(gkdata)) = 0;
        
   
        % reflect the data to sample only 1/2 k-space: symmetry of FT
   %     gkdata(end/2+1: end, :, end/2+1:end) = gkdata( end/2:-1:1, :, end/2+1:end);
   %     gkdata(1:end/2, :, 1:end/2) = gkdata(1:end/2, :, end:-1:end/2+1);
        
        
        subplot(224)
        lightbox(log(abs(gkdata)),[0 7],[]);
        title(sprintf('Log of re-grided k-data, coil %d', c))
        
        
        %  do the FFT :
        gkdata = gkdata ./ (dist+0.001);
        im = fft3d((gkdata ));
        
       %im = xcp;
       
        % combine the coils:
        im_all = abs(im_all) + i*abs(im);
        subplot(223)
        lightbox(abs(im_all)); title('gridding recon')        
        colormap parula

        drawnow
    end
    %%
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
    
    
    
    fprintf('\n ... Done! \n');
    %
    figure(2);
    ov([],abs(im_all), dim3(1)/2, dim3(2)/2, dim3(3)/2, 0);    
    title('Ortho views of gridded recon')

    figure(3)
    ov([],log(abs(gkdata)), dim3(1)/2+1, dim3(2)/2+1, dim3(3)/2+1, 0);
    title('Ortho views of regridded k-space (log)')
    
    
%     figure(3);
%     ov([],abs(xcp), dim3(1)/2, dim3(2)/2, dim3(3)/2, 0);
%     title('Ortho views of NUFFT recon')
    
    
end
%%
return
%%
hdr.dim(1:5) = [4 MYDIM MYDIM nslices ];
xdim = dim(1);
ydim = dim(2);
zdim = dim(3);
tdim = dim(4);
hdr.dim(2:5) = dim;

hdrRaw = ge_pfilehdr(strFileRaw);

hdr.datatype = 4;
hdr.bitpix = 16;

figure; lightbox(abs(imseries)); colormap jet
%{
 figure(2);
 zfov = nslices*slthick/2;
 for n=1:nslices
%       subplot(221)
%      set(gca,'Position',[ 0.05    0.55    0.5  0.5*fov/zfov ]);
%      subplot(222)
%      set(gca,'Position',[ 0.55    0.55    0.5  0.5*fov/zfov ]);
%      subplot(223)
%      set(gca,'Position',[ 0.05    0.05    0.5  0.5 ]);
     ov([],abs(im_all),MYDIM/2,MYDIM/2,n,0);
%    drawnow
     pause(0.5)
 end
%}