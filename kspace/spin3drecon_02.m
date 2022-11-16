% spin3drecon_02
%
% Luis Hernandez-Garcia @UM 2013
%
% 1- This script reads FIDs from the GE scanner acquired with aSPINS
% trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot
% 3 - it re-grids the k-space data into a cartesian grid
% 4 - reconstructs the image by a 3D FFT
%

dt = 4e-6; %seconds
gamma = 267.5 * 1e2 / (2*pi) ; % Hz/Gauss
MYDIM = 64;

myPfile = 'P_ssfptest'

[kdata np nt] = read_raw_3d(myPfile ,  1 );

kdata = reshape(kdata, np, nt);
imagesc(abs(kdata'))

% remove a few points from the begining in case there is a timing offset
Nclip = 1;
kdata = kdata(Nclip+1:end,:);
np = np-Nclip;

k= load('ktraj.txt');
g = load('grad.txt');
rots = load('rotmats.txt');


%kmax = norm(max(k));
% make a filter:
% myfilter = linspace(0,kmax,size(kdata,1))';
%for n=1size(kdata,2)
%    %   kdata(:,n) = kdata(:,n) .* myfilter;
%end

% use the k-space trajectory calculated from the gradients,
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;

% must resample it to match the data acquisition rate
k3 = zeros(np,3);
for n=1:3
    k3(:,n) = resample(k2(:,n), np, size(k2,1));
end

%Now do the reconstruction
% step 1- regrid the k-space data:

Kmax = max(abs(k2(:)));
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM));
ctr = 1;
kxs = []; kys=[]; kzs=[];
for(n=1:3:length(rots))
    m = rots(n:n+2,:);
    kr = m*k3';
    %
    subplot(221)
    plot3(kr(1,:), kr(2,:), kr(3,:));
    hold on
    plot3(krp(1,:), krp(2,:), krp(3,:),'r'); hold off
    axis([-1 1 -1 1 -1 1]*5)
    subplot(222)
    plot(abs(kdata(:,ctr)));
    
    drawnow
    %
    kxs = [kxs kr(1,:)'];
    kys = [kys kr(2,:)'];
    kzs = [kzs kr(3,:)'];
    ctr = ctr +1;
end

% bad data here:
badlist =  cumsum(16*ones(14,1) );
badlist = sort([badlist; badlist-1]);

kdata(:,badlist) = [];
kxs(:,badlist) = [];
kys(:,badlist) = [];
kzs(:,badlist) = [];

kxs = kxs(:)';
kys = kys(:)';
kzs = kzs(:)';

gkdata = zeros(MYDIM, MYDIM, MYDIM);
% gkdata = griddata(kxs, kys, kzs, kdata(:), kx, ky, kz);
%gkdata = interp3(kxs, kys, kzs, kdata(:), kx, ky, kz);
tmpk = kdata(:)';
F = TriScatteredInterp(kxs', kys', kzs', tmpk');
gkdata = F(kx,ky,kz);

%gkdata = gkdata(1:2:end, 1:2:end, 1:2:end);

subplot(223)
lightbox(abs(gkdata),[0 200],[]); title('The re-gridded k-space image')

% step 2: do the FFT recon:
% where the interpolation returned NaN, we must put in zeros
gkdata(isnan(gkdata)) = 0;
im = fftshift(ifftn(fftshift(gkdata)));

%im=im(end/4:3*end/4, end/4:3*end/4, end/4:3*end/4);
%im=im(end/4:3*end/4, end/4:3*end/4, end/4:3*end/4);

subplot(224)
lightbox(abs(im)); title('The reconned image')

fprintf('\n ... Done! \n');
