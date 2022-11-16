FOV = 0.2 ; % m
MYDIM=40;
xres = FOV/MYDIM; % m
gamma = 267.5 * 1e2 % rad/s/Gauss
Tmax = 0.005;  %s
dt = 1/(125e3) ; % s
Kmax = 1/(2*xres);
deltaK = 1/FOV;        % resolution of k-space (separation between turns)
Ntheta = 2*Kmax/deltaK;  % N. turns of the in in-plane spiral required for critical sampling in radial direction
Nphi =   Ntheta;

Nz = 5;
Ntheta = 15;

t=linspace(0,Tmax, Tmax/dt);
Npts = length(t);


THETA = linspace(0, Ntheta*2*pi, Npts);
%PHI = linspace(-Ntheta*pi/2, Ntheta*pi/2, Npts);
R = linspace(0,Kmax^2,Npts);
R = (R).^(1/2);

kxs = []; kys = [];kzs=[];
figure(1)
%subplot(221)
%for dphi = PHI
 for z = linspace(-Kmax,Kmax,Nz)
    
    % [kx, ky, kz] = sph2cart(THETA , dphi*ones(size(R)) , R);
    % [kx, ky, kz] = sph2cart(THETA ,sin(PHI) , R);
    % [kx, ky, kz] = sph2cart(THETA ,(PHI) +dphi , R);
    kx = real(R.*exp(-i*THETA));
    ky = imag(R.*exp(-i*THETA));
    kz = ones(size(kx))*z;
    
    kxs = [kxs; kx'];
    kys = [kys; ky'];
    kzs = [kzs; kz'];
    
    plot3(kx,ky,kz);%axis([-100 100 -100 100 -100 100]); axis square
    drawnow
    pause(0.1)
    hold on
    
    
end
hold off
return
% plot gradients for last shot
Gx = diff(kx) / dt / gamma;
Gy = diff(ky) / dt / gamma;
Gz = diff(kz) / dt / gamma;

xslew = diff(Gx)/dt; % this should come out in gauss/m/s
yslew = diff(Gy)/dt; %
zslew = diff(Gz)/dt; %

subplot(222)
plot(Gx,'r'); hold on
plot(Gy,'g');
plot(Gz,'b');hold off

figure(2)
%% Now simulate an object

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

signal = interp3(kx, ky, kz, kim, kxs, kys, kzs);
figure(1)
subplot(212)
plot(real(signal)); hold on;
plot(imag(signal),'r');
title('All the echoes concatenated')

%Now do the reconstruction
% step 1- regrid the k-space data:

kim2 = griddata3(kxs, kys, kzs, signal, kx, ky, kz);
% kim2 = interp3(kxs, kys, kzs, signal, kx, ky, kz);

figure(2); subplot(224)
lightbox(abs(kim2),[0 1000],[]); title('The re-gridded k-space image')

% step 2: do the FFT recon:
% where the interpolation returned NaN, we must put in zeros
kim2(isnan(kim2)) = 0;
im2 = fftshift(ifftn(fftshift(kim2)));

subplot(223)
lightbox(abs(im2)); title('The reconned image')
