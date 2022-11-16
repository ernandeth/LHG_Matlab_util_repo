FOV = 0.24 ; % m
MYDIM = 64;
THETA_accel = 0.85;
R_accel = 0.75;

xres = FOV/MYDIM; % m
gamma_rad = 267.5 * 1e2     % rad/s/Gauss
gamma_hz = gamma_rad/2/pi;  % Hz/Gauss
Tmax = 0.038;  %s
dt = 4e-6 ; % s
Kmax = 1/(xres);
deltaK = 1/FOV;        % resolution of k-space (separation between turns)
Ntheta = Kmax/deltaK;  % N. turns of the in in-plane spiral required for critical sampling in radial direction
Ntheta = THETA_accel*Ntheta/2;     % since we are going spiral in and spiral out, we need half as many turns in each direction
Nphi =   Ntheta/2;
%Nphi = 0;
Nramp = 800;


t=linspace(0,Tmax, Tmax/dt);
Npts = length(t);


%THETA = linspace(-Ntheta*pi, Ntheta*pi, Npts);
%sfactor = 0.75*pi/2;
%THETA = atan(linspace(-sfactor, sfactor, Npts));
sfactor = pi/3;
THETA = tanh(linspace(-sfactor, sfactor, Npts));
THETA = abs(THETA) * Ntheta*pi/max(abs(THETA));

THETA(end/2 +1:end) = pi + THETA(end/2:-1:1);


PHI = linspace(-Nphi*pi, Nphi*pi, Npts);
R = (abs(linspace(-1, 1, Npts)).^R_accel) ;
R = Kmax*(R);

% smooth out beginning by ramping out to Kmax
R(1:Nramp) = (tanh(linspace(0,2*pi,Nramp))) .* R(1:Nramp);

% reflect the spiral in to get the spiral out
R(end:-1:end/2+1) = R(1:end/2);
R(end) = 0; %-R(end-2);

R = R*Kmax/max(R);  % 1/m


kxs = []; kys = [];kzs=[];
figure(1)
    
[kx, ky, kz] = sph2cart(THETA , zeros(size(R)) , R);
kx = kx'; ky = ky'; kz=kz';

% % fix discontinuity near the center of k-space
Nbuf = 20;
p1 = kx(end/2-Nbuf);
p2 = kx(end/2+Nbuf);
kx(end/2-Nbuf:end/2+Nbuf) = linspace(p1,p2,2*Nbuf+1);

Nbuf = 20;
p1 = ky(end/2-Nbuf);
p2 = ky(end/2+Nbuf);
ky(end/2-Nbuf:end/2+Nbuf) = linspace(p1,p2,2*Nbuf+1);

kx(end/2-200:end/2+200) = smooth(kx(end/2-200:end/2+200),10);
ky(end/2-200:end/2+200) = smooth(ky(end/2-200:end/2+200),10);

kx(end/2-200:end/2+200) = smooth(kx(end/2-200:end/2+200),10);
ky(end/2-200:end/2+200) = smooth(ky(end/2-200:end/2+200),10);

%subplot(221)
cla
plot(kx(1:end/2),ky(1:end/2),'b');
hold on
plot(kx(end/2+1:end),ky(end/2+1:end),'r');
hold off
   
drawnow

% use this code to make the "nautilus spiral"
%{
dphi = 1.5*pi*(linspace(0,1,length(kx)))
for n=1:length(kx)
    rotmat = [
        sin(dphi(n)) 0 cos(dphi(n));
        0   1     0;
        -cos(dphi(n)) 0 sin(dphi(n));
        ];
    
    K = [kx(n) ky(n) kz(n)]* rotmat;
    kx(n) = K(1); ky(n) = K(2); kz(n) = K(3);
end
cla


plot3(kx(1:end/2),ky(1:end/2),kz(1:end/2),'b');
axis([-1 1 -1 1 -1 1]*400); axis square
hold on
plot3(kx(end/2+1:end),ky(end/2+1:end),kz(end/2+1:end),'r');
hold off
%}

   
dphi=[];
%for dphi = deg2rad(117)*[0:Nphi]
golden_angle = deg2rad(137.5);
%for dphi = pi*(linspace(0,1,Nphi))
for dphi =  golden_angle*[0:Nphi-1]
     
    rotmat = [
        cos(dphi) 0 sin(dphi);
        0   1     0;
        -sin(dphi) 0 cos(dphi);
        ];
    
    rotmat2 = [
        1   0   0;
        0   cos(dphi)   sin(dphi);
        0   -sin(dphi)  cos(dphi) ;
        ];
    
    
%    [kx, ky, kz] = sph2cart(THETA , dphi*zeros(size(R)) , R);
%    kx = kx'; ky = ky'; kz=kz';

   
    K = [kx ky kz]* (rotmat*rotmat2);
    
    %kx = K(:,1); ky = K(:,2); kz = K(:,3);
    kxs = [kxs; K(:,1)];
    kys = [kys; K(:,2)];
    kzs = [kzs; K(:,3)];
    
    %
    figure(1)
    plot3(K(1:end/2,1),K(1:end/2,2),K(1:end/2, 3),'b');
    hold on
    plot3(K(end/2+1:end, 1), K(end/2+1:end,2),K(end/2+1:end,3),'r');
    line([-400 400], [0 0] , [0 0])
    line( [0 0] , [-400 400],[0 0])
    line( [0 0] , [0 0], [-400 400])
    hold off
    axis([-1 1 -1 1 -1 1]*400); axis square
    

drawnow
    %}
    
    
end
hold off
figure(2)
subplot(311)
plot(diff(kxs)); 
subplot(312)
plot(diff(kys));
subplot(313)
plot(diff(kzs));
hold off

%% plot gradients for last shot
Gx = diff(kx) / dt / gamma_hz;  % (1/m) / s / Hz/Gauss = Gauss/m
Gy = diff(ky) / dt / gamma_hz;
Gz = diff(kz) / dt / gamma_hz;



% Target slew rates are about 
% 150 - 200 mT/cm/ms
%  (x  1e5) if  mT/m/s
xslew = diff(Gx)/dt*1e-4; % this should come out in mT/m/ms
yslew = diff(Gy)/dt*1e-4; %
zslew = diff(Gz)/dt*1e-4; %

subplot(222)
plot(Gx,'r'); hold on
plot(Gy,'g');
plot(Gz,'b');hold off
title('Gradients (G/m)')

subplot(223)
plot(xslew,'r'); hold on
plot(yslew,'g');
plot(zslew,'b');hold off
title('Slew Rates (mT/m/ms)')

% some metrics about k-space sampling
dRdTH = diff(R)./diff(THETA);
subplot(224)
plot(dRdTH*2*pi); title('2\pi x dR/d\theta')

%return

%% Now simulate an object
figure(2)
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
title('All the echoes overlaid')

%Now do the reconstruction
% step 1- regrid the k-space data:

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
