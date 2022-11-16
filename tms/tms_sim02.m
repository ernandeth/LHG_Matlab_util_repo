function [Ex, Ey, Ez] = tms_sim02(wire, current, sigma, FOV);
% function tms_sim02(wire, current, sigma);

plot3(wire(:,1), wire(:,2), wire(:,3),'r')
hold on , plot3(0,0,0,'ro')
axis([-1 1 -1 1 -1 1])

% Allocate space and define some conductivity space
Nvox = size(sigma,1);
zview = round(Nvox/2);

dx = FOV/Nvox;

sxn = sigma;
syn = sigma;
szn = sigma;

sxd = sigma;
syd = sigma;
szd = sigma;


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

tic
% Calculate E field from coil alone:
[E Ex Ey Ez] =  Efield(current, wire, floor([Nvox Nvox Nvox]/2));
% save Ea.mat E Ex Ey Ez
%load Ea.mat
% figure
% imagesc(squeeze(Ex(:,zview,:))); colorbar
% title ('E_A_y only')
toc
% now get the scalar potential contribution to the E field

% Step 1:  get boundary conditions of Phi 
% based on s1*E1*n = -s2*E2*n at boundaries
% --->  dPhi/dNormal = (s1*E1 + s2*E2) / (s1+s2)
tic
% numerator (s1- s2)* E2
sxn(1:end-1, :, :) = sigma(1:end-1, :, :) - sigma(2:end, :, :) ;
syn(:, 1:end-1, :) = sigma(:, 1:end-1, :) - sigma(:, 2:end, :) ;
szn(:, :, 1:end-1) = sigma(:, :, 1:end-1) - sigma(:, :, 2:end) ;

% denominator: s1 + s2
sxd(1:end-1, :, :) = sigma(1:end-1, :, :) + sigma(2:end, :, :) ;
syd(:, 1:end-1, :) = sigma(:, 1:end-1, :) + sigma(:, 2:end, :) ;
szd(:, :, 1:end-1) = sigma(:, :, 1:end-1) + sigma(:, :, 2:end) ;

% dPhi/dNormal
dPhi_dx = -Ex .* sxn ./ sxd;
dPhi_dy = -Ey .* syn ./ syd;
dPhi_dz = -Ez .* szn ./ szd;

dPhi_dx(find(isnan(dPhi_dx))) = 0;
dPhi_dy(find(isnan(dPhi_dy))) = 0;
dPhi_dz(find(isnan(dPhi_dz))) = 0;

dPhi_dx(find(isinf(dPhi_dx))) = 0;
dPhi_dy(find(isinf(dPhi_dy))) = 0;
dPhi_dz(find(isinf(dPhi_dz))) = 0;

% integrate dPhi/dNormal to get a potential (Phi) map for initial condition
% Phi = cumsum(dPhi_dx,1) + cumsum(dPhi_dy,2) + cumsum(dPhi_dz,3);
Phix = cumsum(dPhi_dx,1) ;
Phiy = cumsum(dPhi_dy,2) ;
Phiz = cumsum(dPhi_dz,3) ;
Phi = sqrt ( Phix.^2 + Phiy.^2 + Phiz.^2);

% Phi = sigma;
toc
% Step 2: Now solve Laplace equation numerically given inital values of
% Phi map.
old = Phi;
new = old;
new(:,:,:)=0;

energy = sum(sum(sum(Phi)));
tic
tol = 0.05*energy;
MSE = 10;
% figure
while (MSE >= tol)
    new(2:Nvox-1, 2:Nvox-1 , 2:Nvox-1) = ...
       (old(1:Nvox-2, 2:Nvox-1, 2:Nvox-1) + old(3:Nvox, 2:Nvox-1, 2:Nvox-1) + ...
        old(2:Nvox-1, 1:Nvox-2, 2:Nvox-1) + old(2:Nvox-1, 3:Nvox, 2:Nvox-1) + ...
        old(2:Nvox-1, 2:Nvox-1, 1:Nvox-2) + old(2:Nvox-1, 2:Nvox-1, 3:Nvox)) / 6;

    % Keep the surface potential constant, despite the iterations.
    new(surfVoxels) = Phi(surfVoxels);

    MSE = abs(old - new);,
    MSE = sum(reshape(MSE, Nvox^3,1));
    old = new;
    fprintf('\rMSE = %f', MSE/energy);
%     imagesc(new(:,:,zview));
%     drawnow
end
toc
Phi = new;
% Step 3: E_phi is given by Maxwell's equation E = -grad(Phi)
[Ephix Ephiy Ephiz] = gradient(Phi,dx,dx,dx);

% Now add both contributions to the Electric field.
Ex_total = Ex - Ephix;
Ey_total = Ey - Ephiy;
Ez_total = Ez - Ephiz;
% 
% subplot(211)
% quiver(Ex(:,:,zview), Ey(:,:,zview));
% subplot(212)
% quiver(Ex_total(:,:,zview), Ey_total(:,:,zview));
% %imagesc(squeeze(Ex_total(:,:,zview))); colorbar
% title(sprintf('slice %d', zview));
% 
% imagesc(squeeze(Ex_total(:,:,zview))); colorbar
% imagesc(squeeze(Ez_total(:,:,zview))); colorbar
% imagesc(squeeze(Phi(:,:,zview))); colorbar

return