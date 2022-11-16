dbstop if error
close all
deg=0 
current =1;
config='hemisphere'
r = 0.4
R = 0.8
wire=[];
axis xy

FOV = 0.24; % m

% Make the coils:
loopSegs = linspace(0,2*pi-2*pi/20,20);
wx = r*sin(loopSegs);
wy = -r*cos(loopSegs);
wz = R*ones(size(wx));
wire =[wx' wy' wz';]
% loop =[wx' wy' wz';];
% %loop =[wx' -0.4-wy' wz';];
% %     wx' 0.4+wy' wz'];
% 
% % make more loops by rotation about the y axis
% % for theta= 0:pi/2:pi/2
% for theta= 0:0
%     rotmat = ...
%         [cos(theta)   0         -sin(theta) ;
%         0             1         0 ;
%         sin(theta)   0         cos(theta)];
%     loop = loop * rotmat;
%     wire = [wire ; loop];
% end

% for debugging, use just a line:
% wire = [-0.5 0 0.8;
%     0.5 0 0.8];

plot3(wire(:,1), wire(:,2), wire(:,3),'r')
%hold on,plot(wire(:,1), wire(:,2),'r')
hold on , plot3(0,0,0,'ro')
axis([-1 1 -1 1 -1 1])

% Allocate space and define some conductivity space
Nvox = 64;
zview = 30;
sigma = zeros(Nvox,Nvox,Nvox);

dx = FOV/(Nvox-1);

sxn = sigma;
syn = sigma;
szn = sigma;

sxd = sigma;
syd = sigma;
szd = sigma;

sigma(22:43, 22:43, 22:43) = 0.7; % Siemens/meter or kg^?1ám^?3ás^3áA^2
%sigma(27:37, 27:37, 27:37) = 0; 

%load sigma

[dsdx, dsdy, dsdz] = gradient(sigma);
surfVoxels = find( dsdx | dsdy | dsdz);

Ex_total = zeros(Nvox,Nvox,Nvox);
Ey_total = zeros(Nvox,Nvox,Nvox);
Ez_total = zeros(Nvox,Nvox,Nvox);

dPhi_dx = zeros(Nvox,Nvox,Nvox);
dPhi_dy = zeros(Nvox,Nvox,Nvox);
dPhi_dz = zeros(Nvox,Nvox,Nvox);

Phix = zeros(Nvox,Nvox,Nvox);
Phiy = zeros(Nvox,Nvox,Nvox);
Phiz = zeros(Nvox,Nvox,Nvox);

% Calculate E field from coil alone:
[E Ex Ey Ez] =  Efield(current, wire, floor([Nvox Nvox Nvox]/2));
save Ea.mat E Ex Ey Ez
%load Ea.mat
%figure
% for zview = 1:64
%     zview
%     subplot 311
%     quiver(Ex(1:3:end,1:3:end,zview), Ey(1:3:end,1:3:end,zview));
%     subplot 312
%     imagesc(Ex(:,:,zview));
%     subplot 313
%     imagesc(Ey(:,:,zview));
%     pause
% end

% for zview =1:Nvox;
%     quiver(Ex(:,:,zview), Ey(:,:,zview));
%     pause
% end

% now get the scalar potential contribution to the E field

% Step 1:  get boundary conditions of Phi 
% based on s1*E1*n = -s2*E2*n at boundaries
% --->  dPhi/dNormal = (s1*E1 + s2*E2) / (s1+s2)

% numerator (s1- s2)* E2
sxn(1:end-1, :, :) = -sigma(1:end-1, :, :) + sigma(2:end, :, :) ;
syn(:, 1:end-1, :) = -sigma(:, 1:end-1, :) + sigma(:, 2:end, :) ;
szn(:, :, 1:end-1) = -sigma(:, :, 1:end-1) + sigma(:, :, 2:end) ;

% denominator: s1 + s2
sxd(1:end-1, :, :) = sigma(1:end-1, :, :) + sigma(2:end, :, :) ;
syd(:, 1:end-1, :) = sigma(:, 1:end-1, :) + sigma(:, 2:end, :) ;
szd(:, :, 1:end-1) = sigma(:, :, 1:end-1) + sigma(:, :, 2:end) ;

% dPhi/dNormal
dPhi_dx = abs(Ex) .* sxn ./ sxd;
dPhi_dy = abs(Ey) .* syn ./ syd;
dPhi_dz = abs(Ez) .* szn ./ szd;

dPhi_dx(find(isnan(dPhi_dx))) = 0;
dPhi_dy(find(isnan(dPhi_dy))) = 0;
dPhi_dz(find(isnan(dPhi_dz))) = 0;

dPhi_dx(find(isinf(dPhi_dx))) = 0;
dPhi_dy(find(isinf(dPhi_dy))) = 0;
dPhi_dz(find(isinf(dPhi_dz))) = 0;

% integrate dPhi/dNormal to get a potential (Phi) map for initial condition
Phix = cumsum(dPhi_dx,1) *dx ;
Phiy = cumsum(dPhi_dy,2) *dx;
Phiz = cumsum(dPhi_dz,3) *dx;
% Phix (2:end, :, :) = Phix(1:end-1,: , :) + dPhi_dx(1:end-1,: , :) *dx;
% Phiy (:, 2:end, :) = Phiy(:, 1:end-1, :) + dPhi_dy(:, 1:end-1, :) *dx;
% Phiz (:, :, 2:end) = Phiz(:, :, 1:end-1) + dPhi_dz(:, :, 1:end-1) *dx;

%Phi = sqrt ( Phix.^2 + Phiy.^2 + Phiz.^2);
Phi = (Phix + Phiy + Phiz)/3;
Phi(:,:,end) = 0;
Phi(:,end,:) = 0;
Phi(end,:,:) = 0;


% Phi = sigma;

% Step 2: Now solve Laplace equation numerically given inital values of
% Phi map.
old = Phi;
new = old;
new(:,:,:)=0;

energy = sum(Phi(:).^2);
colormap gray
tol = 0.000001*energy;
MSE = 10;
figure

while (MSE >= tol)
    new(2:Nvox-1, 2:Nvox-1 , 2:Nvox-1) = ...
       (old(1:Nvox-2, 2:Nvox-1, 2:Nvox-1) + old(3:Nvox, 2:Nvox-1, 2:Nvox-1) + ...
        old(2:Nvox-1, 1:Nvox-2, 2:Nvox-1) + old(2:Nvox-1, 3:Nvox, 2:Nvox-1) + ...
        old(2:Nvox-1, 2:Nvox-1, 1:Nvox-2) + old(2:Nvox-1, 2:Nvox-1, 3:Nvox)) / 6;

    % Keep the surface potential constant, despite the iterations.
    new(surfVoxels) = Phi(surfVoxels);

    MSE = abs(old - new);,
    MSE = sum(MSE(:).^2);
    old = new;
    imagesc(squeeze(new(:,zview,:)));
    fprintf('\rMSE = %f', MSE/energy);
    axis xy; title('sagittal Phi Map updating')
    drawnow
end

Phi = new;
% Step 3: E_phi is given by Maxwell's equation E = -grad(Phi)
[Ephix Ephiy Ephiz] = gradient(Phi,dx,dx,dx);

% Now add both contributions to the Electric field.
Ex_total = Ex - Ephix;
Ey_total = Ey - Ephiy;
Ez_total = Ez - Ephiz;

figure
subplot(211)
quiver(Ex(1:3:end,1:3:end,zview), Ey(1:3:end,1:3:end,zview));
axis xy; title(sprintf('E_A slice %d transverse', zview));
subplot(212)
quiver(Ex_total(1:3:end,1:3:end,zview), Ey_total(1:3:end,1:3:end,zview));
%imagesc(squeeze(Ex_total(:,:,zview))); colorbar
axis xy; title(sprintf('Etotal slice %d transverse', zview));

figure
subplot 311
imagesc(squeeze(Ex_total(:,:,zview))); colorbar
axis xy; title('Ex transverse')
subplot 312
imagesc(squeeze(Ez_total(:,zview,:))); colorbar
axis xy; title('Ez sagittal')
subplot 313
imagesc(squeeze(Phi(:,:,zview))); colorbar
axis xy; title('Phi tranxverse')

figure
subplot 311, plot(squeeze(Ex_total(zview, zview,:)));
axis xy; title('Ex z-axis')
subplot 312, plot(squeeze(Ey_total(zview, zview,:)));
axis xy; title('Ey z-axis')
subplot 313, plot(squeeze(Phi(zview, zview,:)));
axis xy; title('Phi z-axis')

figure
subplot 311, plot(squeeze(Ex(:,zview, zview)));
axis xy; title('E_ax x-axis')
subplot 312, plot(squeeze(Ey(:,zview, zview)));
axis xy; title('E_ay x-axis')
subplot 313, plot(squeeze(Phi(:,zview, zview)));
axis xy; title('Phi x-axis')

figure
subplot 311, plot(squeeze(Ex_total(:,zview, zview)));
axis xy; title('Ex x-axis')
subplot 312, plot(squeeze(Ey_total(:,zview, zview)));
axis xy; title('Ey x-axis')
subplot 313, plot(squeeze(Phi(:,zview, zview)));
axis xy; title('Phi x-axis')

figure
subplot 311
imagesc(squeeze(dPhi_dx(:,:,zview))); colorbar
axis xy; title('dPhi/dx transverse')
subplot 312
imagesc(squeeze(dPhi_dy(:,:,zview))); colorbar
axis xy; title('dPhi/dy transverse')
subplot 313
imagesc(squeeze(Phi(:,:,zview))); colorbar
axis xy; title('Phi tranxverse')

figure
subplot 311
imagesc(squeeze(Phix(:,:,zview))); colorbar
axis xy; title('Phix transverse')
subplot 312
imagesc(squeeze(Phiy(:,:,zview))); colorbar
axis xy; title('Phiy transverse')
subplot 313
imagesc(squeeze(Phi(:,:,zview))); colorbar
axis xy; title('Phi tranxverse')

figure
subplot 311
plot(squeeze(dPhi_dx(:,zview,zview))); colorbar
axis xy; title('dPhi/dx x axis')
subplot 312
plot(squeeze(dPhi_dy(zview,:,zview))); colorbar
axis xy; title('dPhi/dy y axis')
subplot 313
plot(squeeze(Phi(:,zview,zview))); colorbar
axis xy; title('Phi x axis')

figure
subplot 311
plot(squeeze(Phix(:,zview,zview))); colorbar
axis xy; title('Phix x axis')
subplot 312
plot(squeeze(Phiy(zview,:,zview))); colorbar
axis xy; title('Phiy y axis')
subplot 313
plot(squeeze(Phi(zview,:,zview))); colorbar
axis xy; title('Phi tranxverse')
