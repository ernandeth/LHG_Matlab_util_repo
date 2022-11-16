FOV = 24 ; % cm
MYDIM=64;
xres = FOV/MYDIM; % m
gamma = 267.5 * 1e2 / (2*pi) ; % Hz/Gauss
Tmax = 0.01;  %s
dt = 4e-6 ; % s
dt = 8e-6;
Kmax = 1/(xres);
deltaK = 1/(2*FOV);        % resolution of k-space (separation between turns)

Nshots = 25;

% new 10.2
Ntheta = Kmax/deltaK;
Nphi = Kmax/deltaK;


t=linspace(0,Tmax, Tmax/dt);
Npts = length(t);


THETA = linspace(0, 2*Ntheta*pi/Nshots, Npts);
PHI =   linspace(0, 2*Nphi*pi/Nshots, Npts);
R = linspace(0, Kmax^(1), Npts);
R = (R).^(1);


dphi = 2*pi/Nshots;
dtheta = 2*pi/(Nshots);



figure(1)
subplot(221)

option=2
if option==1
    kxs = []; kys = [];kzs=[];
    for n=0:(Nshots)-1
            
        PHI = PHI+ dphi;
            
        for m=0:(Nshots)-1
            
            THETA = THETA + dtheta;
            
            [kx, ky, kz] = sph2cart(PHI, THETA, R);
            
            kxs = [kxs; kx'];
            kys = [kys; ky'];
            kzs = [kzs; kz'];
            
            %plot3(kx,ky,kz);axis([-2 2 -2 2 -2 2]); axis square
            %drawnow
            %pause(0.05)
            %hold on
        end
        
        
    end
    hold off
    
else%% alternative: cartesian coordinate version for EPIC prog.
    
    % first spiral:
    [kx, ky, kz] = sph2cart(THETA, PHI, R);
    kxs = []; kys = [];kzs=[];
    plot3(kx,ky,kz , 'g');axis([-2 2 -2 2 -2 2]); axis square
    
    rotz= [
        cos(dtheta)  -sin(dtheta)   0;
        sin(dtheta)  cos(dtheta)    0;
        0    0  1;
        ];
    
    rotx= [
        1   0   0;
        0   cos(dphi)  -sin(dphi);
        0   sin(dphi)  cos(dphi) ;
        ];
    
    R=eye(3);
    for n=0:(Nshots)-1
        
        for m=0:(Nshots)-1
          
            
            
            tmp = R * [kx ; ky ; kz];
            kx1 = tmp(1,:);
            ky1 = tmp(2,:);
            kz1 = tmp(3,:);
            
            kxs = [kxs; kx1'];
            kys = [kys; ky1'];
            kzs = [kzs; kz1'];
            
            %plot3(kx1,ky1,kz1 , 'r');  
            %axis([-2 2 -2 2 -2 2]); axis square
            
            %drawnow
            %pause(0.1)
            %hold on
            R=rotx*R;
        end
        R=rotz*R;
        
    end
    hold off
    
    
    
    
end



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
axis([0 Npts -3 3])
title('Gradient Waveforms')




figure(2)
% Now simulate an object

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

signal = interp3(kx, ky, kz, kim, kxs, kys, kzs,'spline');
figure(1)
subplot(212)
plot(real(signal)); hold on;
plot(imag(signal),'r');
plot(abs(signal),'g');

title('All the echoes concatenated')
hold off

% filter the data
% make a filter:
myfilter = linspace(eps,Kmax,Npts)';
myfilter = repmat(myfilter, Nshots*Nshots, 1);
% signal = signal .* myfilter;

%Now do the reconstruction
% step 1- regrid the k-space data:

kim2 = griddata(kxs, kys, kzs, signal, kx, ky, kz);
% kim2 = interp3(kxs, kys, kzs, signal, kx, ky, kz);

figure(2); subplot(224)
lightbox(abs(kim2),[0 1000],[]); title('The re-gridded k-space image')

% step 2: do the FFT recon:
% where the interpolation returned NaN, we must put in zeros
kim2(isnan(kim2)) = 0;
im2 = fftshift(ifftn(fftshift(kim2)));

subplot(223)
lightbox(abs(im2)); title('The reconned image')

%RMK added for densitiy compenstation
figure(3)
subplot(2,1,1)
plot(abs(kim(:,32,32)))
hold on
plot(abs(kim(32,:,32)),'r')
hold off
subplot(2,1,2)
plot(abs(kim2(:,32,32)))
hold on
plot(abs(kim2(32,:,32)),'r')
hold off

figure(4)
diff1 = (abs(kim(:,32,32))-abs(kim2(:,32,32)))./(abs(kim(:,32,32)))
diff2 = (abs(kim(32,:,32))-abs(kim2(32,:,32)))./(abs(kim(32,:,32)))
plot(diff1);
hold on
plot(diff2,'r');
hold off



%End RMK


fprintf('\n ... Done! \n');
