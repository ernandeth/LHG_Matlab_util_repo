FOV = 0.22 ; % m
MYDIM = 64;
xres = FOV/MYDIM; % m
gamma = 267.5 * 1e2 % rad/s/Gauss
Tmax = 0.020;  %s
dt = 4e-6 ; % s
Kmax = 1/(xres);
deltaK = 1/FOV;        % resolution of k-space (separation between turns)
Ntheta = Kmax/deltaK;  % N. turns of the in in-plane spiral required for critical sampling in radial direction
Ntheta = Ntheta/2;     % since we are going spiral in and spiral out, we need half as many turns in each direction
Nphi =   Ntheta/2;
%Nphi = 0;
Nramp = 200;


t=linspace(0,Tmax, Tmax/dt);
Npts = length(t);


%THETA = linspace(-Ntheta*pi, Ntheta*pi, Npts);
%THETA = tanh(linspace(-sfactor, sfactor, Npts));
sfactor = 0.75*pi/2;
THETA = atan(linspace(-sfactor, sfactor, Npts));
THETA = THETA * Ntheta*pi/max(THETA);
THETA(end/2 +1:end) = pi+THETA(end/2:-1:1);
PHI = linspace(-Nphi*pi, Nphi*pi, Npts);
R = abs(linspace(-1,1,Npts));
R = Kmax*(R);

%% Now the Duyn method:  I can't get the scaling right,
% but we can just scale it by brute force
maxSlew = 15*1e5; % G/m/s
Taq = 2 * (pi^1.5) * MYDIM / 3 / sqrt(maxSlew*xres);  % seconds
a2 = MYDIM * pi / (Taq^0.666667);
a1 = (3/2) * maxSlew / a2;

% scaling with Duyn paper:
%THETA(1:end/2) = a2*(linspace(Taq,0,Npts/2).^(2/3));
%R(1:end/2) = a1*linspace(Taq, 0, Npts/2) .^(1/3);

% That didn't work so well ...? 
% brute force scaling:
THETA(1:end/2) = Ntheta * pi * abs(linspace(1,0,Npts/2).^(2/3));
R(1:end/2) = abs(Kmax* linspace(1,0,Npts/2) .^(1/3));

% spiral in is almost the same as spiral out:
THETA(end/2+1:end) = THETA(end/2:-1:1) - pi;
R(end/2+1:end) = R(end/2:-1:1);

% ramping up nicely:
% R(1:Nramp) = atan(linspace(0,pi/2,Nramp)) .* R(1:Nramp)/atan(pi/2);
% R(end:-1:end-Nramp+1) = R(1:Nramp);
% R = R*Kmax/max(R);


kxs = []; kys = [];kzs=[];
figure(1)

[gx, gy, gz] = sph2cart(THETA , zeros(size(R)) , R);
kx = gamma*cumsum(gx)*dt;
ky = gamma*cumsum(gy)*dt;
kz = gamma*cumsum(gz)*dt;



[kx, ky, kz] = sph2cart(THETA , zeros(size(R)) , R);

cla
% some metrics about k-space sampling
dRdTH = diff(R)./diff(THETA);
subplot(224)
plot(dRdTH*2*pi); title('2\pi x dR/d\theta')

%{
dphi = (linspace(0,2*pi,length(kx)))
for n=1:length(kx)
    rotmat = [
        sin(dphi(n)) 0 cos(dphi(n));
        0   1     0;
        0 -cos(dphi(n)) sin(dphi(n));
        ];
    
    K = [kx(n) ky(n) kz(n)]* rotmat;
    kx(n) = K(1); ky(n) = K(2); kz(n) = K(3);
end

%}
subplot(221)
plot3(kx(1:end/2),ky(1:end/2),kz(1:end/2),'b');
axis([-1 1 -1 1 -1 1]*300); axis square
hold on
plot3(kx(end/2+1:end),ky(end/2+1:end),kz(end/2+1:end),'r');
hold off

K0 = [kx' ky' kz'];
    
for dphi = pi*(linspace(0,1,Nphi))
    
    rotmat = [
        sin(dphi) 0 cos(dphi);
        0   1     0;
        -cos(dphi) 0 sin(dphi);
        ];
    
    
    %[kx, ky, kz] = sph2cart(THETA , dphi*zeros(size(R)) , R);
    
    K = K0 * rotmat;
    
   
    subplot(221)
    plot3(K(1:end/2,1),K(1:end/2,2),K(1:end/2,3),'b');
    hold on
    plot3(K(end/2+1:end,1),K(end/2+1:end,2),K(end/2+1:end,3),'r');
    %hold off
    axis([-1 1 -1 1 -1 1]*350); axis square
    drawnow
    
    
    kxs = [kxs; K(:,1)];
    kys = [kys; K(:,2)];
    kzs = [kzs; K(:,3)];
    
end
hold off

%% plot gradients for last shot
Gx = diff(kx) / dt / gamma;
Gy = diff(ky) / dt / gamma;
Gz = diff(kz) / dt / gamma;


% Target slew rates are about
% 15 - 20 G/cm/ms
%  (x  1e5) if  G/m/s
xslew = diff(Gx)/dt*1e-4; % this should come out in mT/m/ms
yslew = diff(Gy)/dt*1e-4; %
zslew = diff(Gz)/dt*1e-4; %

subplot(222)
plot(Gx,'r'); hold on
plot(Gy,'g');
plot(Gz,'b');hold off
title('Gradients')

subplot(223)
plot(xslew,'r'); hold on
plot(yslew,'g');
plot(zslew,'b');hold off
title('Slew Rates')

%return

%% Now simulate an object
figure(2)
%MYDIM = MYDIM*2;
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
subplot(224)
plot(real(signal)'); hold on;
plot(imag(signal)','r');
title('All the echoes concatenated')

%Now do the reconstruction
% step 1- regrid the k-space data:

% MYDIM = MYDIM/2;
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM));

kim2 = griddata(kxs, kys, kzs, signal, kx, ky, kz, 'natural');
% kim2 = interp3(kxs, kys, kzs, signal, kx, ky, kz);

figure(2); subplot(224)
lightbox(abs(kim2),[0 1000],[]); title('The re-gridded k-space image')

% step 2: do the FFT recon:
% where the interpolation returned NaN, we must put in zeros
kim2(isnan(kim2)) = 0;
im2 = fftshift(ifftn(fftshift(kim2)));

subplot(223)
lightbox(abs(im2)); title('The reconned image')


return
%%
opxres = 64
gfov = 22
oprbw = 125

Kmax = opxres / gfov;  %/* cm^-1*/
deltaK = 1 / gfov;
FID_len = (Kmax*Kmax / deltaK/deltaK/4);
FID_len = FID_len*(1.0 / (3.14259*0.5*0.5))
FID_dur = FID_len / oprbw/1e3 %/* seconds*/
