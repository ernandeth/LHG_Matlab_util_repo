% spiral3drecon
%
% Luis Hernandez-Garcia @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - it re-grids the k-space data into a cartesian grid
% 4 - reconstructs the image by a 3D FFT
%

dt = 4e-6; %seconds
gamma = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma/2/pi ; % Hz/Gauss

MYDIM = 64;
slthick = 0.6;
myPfile = 'P08704.7'

% load the data:
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

% first add all the coils together
tmp = zeros(nframes, nslices*ndat);
for n = 1:ncoils     
     tmp = abs(abs(tmp) + i*abs(raw(:,:,n)));    
end
kdata = tmp'/ncoils;  % now kdata is ndat*nslices x nframes

subplot(221)
plot(abs(kdata(:,1))); 
title('Echo train: RMS over coils');
drawnow

% Now reshape into a 3D matrix 
kdata = reshape(kdata,  ndat, nslices, nframes);

fprintf('\nLoading grads...');
g = load('grad.txt');
glen = size(g,1);




for Ndel= 0
    
% use the k-space trajectory calculated from the gradients,
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;

% delays in gradient/acquisition?
%Ndel = 100;
k2 = [zeros(abs(Ndel),3); k2; zeros(abs(Ndel),3)];
k2 = circshift(k2,Ndel,1);
k2 = k2(abs(Ndel)+1:end-abs(Ndel),:);

% load k-space trajectory in 2D
% k2= load('ktraj_cart.txt');

k2 =  k2';

% must resample it to match the data acquisition rate
ks = zeros(3,ndat);
fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
chunk = glen-ndat;
for n=1:3
    tmp = k2(n,:);
    ftmp = fft(tmp);
    ftmp = [ftmp(1:end/2-chunk/2) ftmp(end/2+chunk/2+1:end)];
    ks(n,:) = real(ifft(ftmp)) * ndat/glen;
end


subplot(222)
plot(ks(1,:), ks(2,:));
title('Trajectory k-space traj. from grad file ');

%Now do the reconstruction
% step 1- regrid the k-space data:
fprintf('\nCreating target cartesian grid ...');

% Target grid
Kxymax = max(abs(k2(:)));
Kxymax = MYDIM/fov;  % Kmax
Kzmax = 1/(slthick)/2;
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kxymax,Kxymax,MYDIM), ...
    linspace(-Kxymax,Kxymax,MYDIM), ...
    linspace(-Kzmax,Kzmax,nslices));

% Actually sampled  grid : stack of spirals
% kx and ky are always the same, but kz increases monotonically
% this was not recorded in the grad file
fprintf('\nCreating whole 3D sampled grid ...');
tmp = repmat(ks,  1, nslices);
kxs = tmp(1,:);
kys = tmp(2,:);
kzs = tmp(3,:);

for n=1:nslices
    beg = ndat*(n-1)+1;
    fin = ndat*n;
    kzz = 2*Kzmax *(n-1)/ nslices - Kzmax;
    kzs(beg:fin) = kzz;
end

%%
gkdata = zeros(MYDIM, MYDIM, nslices);

for n=2 % :nframes
    fprintf('\nRegridding frame %d ...', n);
    tmpk = kdata(:,:,n);
    tmpk = tmpk(:);
    
    %tmpk = circshift(tmpk,-100);
    
    gkdata = griddata(kxs(:), kys(:), kzs(:), tmpk(:), kx, ky, kz);
    %gkdata = interp3(kxs(:), kys(:), kzs(:), tmpk(:), kx(:), ky(:), kz(:), 'spline');
    %gkdata = reshape(gkdata, MYDIM, MYDIM, nslices);
    %F = TriScatteredInterp(kxs(:), kys(:), kzs(:), tmpk(:));
%     F = scatteredInterpolant(kxs(:), kys(:), kzs(:), tmpk(:));
%     gkdata = F(kx,ky,kz);
    
    %gkdata = gkdata(1:2:end, 1:2:end, 1:2:end);
    
    % where the infidserpolation returned NaN, we must put in zeros
    gkdata(isnan(gkdata)) = 0;
    %gkdata = gkdata - mean(gkdata(:));
    
    subplot(223)
    lightbox(log(abs(gkdata)),[],[]); title('Log k-space image (re-grided)')
    
    % step 2: do the FFT recon:
 
    im = fft3d(gkdata);
    %im = fftshift(fft2(squeeze(gkdata(:,:,4))));
    
    subplot(224)
    lightbox((abs(im))); title(sprintf('The reconned image, del=%d', Ndel))
    
end
fprintf('\n ... Done! \n');
drawnow
end
