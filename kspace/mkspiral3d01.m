FOV = 0.24 ; % m
xres = 3e-3; % m
gamma = 267.5 * 1e2 % rad/s/Gauss
Nphi = 16
Tmax = 0.020;  %s
dt = 4e-6 ; % s
Kmax = 1/xres;
deltaK = 1/FOV;        % resolution of k-space (separation between turns)
Ntheta = Kmax/deltaK;  % N. turns of the in in-plane spiral 

t=linspace(0,Tmax, Tmax/dt);
Npts = length(t);

THETA = linspace(0, Ntheta*2*pi, Npts);
PHI = linspace(0, pi, Nphi);
R = linspace(0,Kmax.^2,Npts);

dphi2 =linspace(0,2*pi,Npts) ;
clf
kxs=[];kys=[];kzs=[];

figure(1)

% Now simulate an object
MYDIM=40;
im = phantom3d(MYDIM);
kim = fftshift(fftn(fftshift(im)));
subplot(221)
lightbox(im); title('The object')
subplot(222)
lightbox(abs(kim),[0 1000],[]); title('THe K-space version of the object')

% a k-space grid
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM));

figure(2)
for dphi=PHI
    
    x1 = sqrt(R).*sin(THETA)  ;
    y1 = sqrt(R).*cos(THETA) ;
    z1 = 0;
    
    % rotate by dphi along the y axis
    x = x1.*cos(dphi+dphi2) + z1.*sin(dphi+dphi2);
    y = y1;
    z = -x1.*sin(dphi+dphi2) + z1.*cos(dphi+dphi2);
    
    
    Gx = diff(x) / dt / gamma;
    Gy = diff(y) / dt / gamma;
    Gz = diff(z) / dt / gamma;
    
    subplot(131)
  
    plot3(x,y,z)
    
    axis([-Kmax Kmax -Kmax Kmax -Kmax Kmax]);
    axis square
    title ('kspace trajectory')
    drawnow
    
    subplot(132)
    plot(Gx, 'r'); hold on
    plot(Gy,'g');
    plot(Gz,'b'); hold off
    title('gradient waveforms in G/m')
    drawnow
    
    kxs = [kxs x];
    kys = [kys y];
    kzs= [kzs z];
    
end
hold off

xslew = diff(Gx)/dt; % this should come out in gauss/m/s 
yslew = diff(Gy)/dt; %
zslew = diff(Gz)/dt; % 


signal = interp3(kx, ky, kz, kim, kxs, kys, kzs);
subplot(133)
plot(real(signal)); hold on;
plot(imag(signal),'r');
title('All the echoes concatenated')

%Now do the reconstruction
% step 1- regrid the k-space data:

kim2 = griddata3(kxs, kys, kzs, signal, kx, ky, kz);
figure(1)
subplot(224)
lightbox(abs(kim2)); title('The re-gridded k-space image')

% step 2: do the FFT recon:
% where the interpolation returned NaN, we must put in zeros
kim2(isnan(kim2)) = 0;
im2 = fftshift(ifftn(fftshift(kim2))); 

subplot(223)
lightbox(abs(im2)); title('The reconned image')
