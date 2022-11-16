% spin3drecon_01
%
% Luis Hernandez-Garcia @UM 2013
%
% 1- This script reads FIDs from the varian scanner acquired with aSPINS
% trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - it re-grids the k-space data into a cartesian grid
% 4 - reconstructs the image by a 3D FFT
%
clear all;
close all;

dt = 4e-6; %seconds
gamma = 267.5 * 1e2 / (2*pi) ; % Hz/Gauss
MYDIM = 64;

[kdata, np, nt] = ReadVarian2D('fid'); plot(abs(kdata))
kdata = reshape(kdata, np/2, nt);
imagesc(abs(kdata'))

% remove a few points from the begining in case there is a timing offset
Nclip = 5;
kdata = kdata(Nclip+1:end,:);
np = np-Nclip*2;

k= load('ktraj.txt');
g = load('grad.txt');
rots = load('rotmats.txt');

% kmax = norm(max(k));
% % make a filter:
% myfilter = linspace(0,kmax,size(kdata,1))';
% for n=1size(kdata,2)
%     kdata(:,n) = kdata(:,n) .* myfilter;
% end

% use the k-space trajectory calculated from the gradients, 
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;
k2 = k;

% must resample the kspace trajectory to match the data acquisition rate
k3 = zeros(np/2,3);
for n=1:3
    k3(:,n) = resample(k2(:,n), np/2, size(k2,1));
end

%Now do the reconstruction
% step 1- regrid the k-space data:

Kmax = max(abs(k3(:)));
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM));

ctr = 1;
kxs = []; kys=[]; kzs=[];
for(n=1:3:length(rots))
    
    m = rots(n:n+2,:);
    kr = m*k3';
    krp = [[0 0 0]' kr(:,end)];
    gr = m*g';
    %
    subplot(221)
    plot3(kr(1,:), kr(2,:), kr(3,:));
    hold on
    plot3(krp(1,:), krp(2,:), krp(3,:),'r'); hold off
    
    axis([-1 1 -1 1 -1 1]*5)
    subplot(222)
    plot(abs(kdata(:,ctr)));


    subplot(223)
    plot(gr')
    %
    pause(0.1)
    drawnow
    
    kxs = [kxs kr(1,:)];
    kys = [kys kr(2,:)];
    kzs = [kzs kr(3,:)];
    
    ctr = ctr +1;
end

% return
gkdata = zeros(MYDIM, MYDIM, MYDIM);
%gkdata = griddata(kxs, kys, kzs, kdata(:), kx, ky, kz);
%gkdata = griddata(kxs(1:end/2), kys(1:end/2), kzs(1:end/2), kdata((1:end/2)), kx, ky, kz);

% gkdata = interp3(kxs, kys, kzs, kdata(:)', kx, ky, kz, 'spline');
% gkdata = gkdata(1:2:end, 1:2:end, 1:2:end);
% tmp = scatteredInterpolant(kxs, kys, kxs, kdata(:))'
% gkdata = tmp(kx,ky,kz);

tmpk = kdata(:)';
F = TriScatteredInterp(kxs', kys', kzs', tmpk')
gkdata = F(kx,ky,kz);

subplot(223)
lightbox(abs(gkdata),[0 100],[]); title('The re-gridded k-space image')

% step 2: do the FFT recon:

% set to zero those locations where the interpolation returned NaN,
gkdata(isnan(gkdata)) = 0;

im = fftshift(ifftn(fftshift(gkdata)));

%im=im(end/4:3*end/4, end/4:3*end/4, end/4:3*end/4); 
%im=im(end/4:3*end/4, end/4:3*end/4, end/4:3*end/4); 

%subplot(224)
figure
%lightbox(abs(im(:,:,25:39))); title('The reconned image')
lightbox(abs(im)); title('The reconned image')

fprintf('\n ... Done! \n');
