function [signal b1map ] = sim_kspace3d_sense(kxs, kys, kzs)
% function [signal b1maps ] = sim_kspace3d_sense(kxs, kys, kzs)
%
% (c) Luis Hernandez-Garcia @UM 2022
%
% k space signal simulator.
%
% Generates a 3d phantom and Ncoils arbtray sensitivity maps
% then simulates the signal observed by the specified k-space trajectory
% Neglecting decay , off-resonance ... or anything else
% Returns the observed signal at each of the Ncoils coils
% along with the sensitivity maps
%

%%
MYDIM = 65;
Ncoils = 8;
ks = [kxs, kys, kzs];
Kmax = max(abs(ks(:)));

% Make some very simplistic complex sensitivity maps
b1map = zeros(MYDIM, MYDIM, MYDIM, Ncoils);
kim = b1map;
fieldw = 600;
for m=1:MYDIM
    for n=1:MYDIM
        for p=1:MYDIM
            b1map(m,n,p,1) = 0.15 + exp(-i*pi/9 ) * exp(-((m-35).^2 + (n-60).^2 + (p-35).^2 )/ fieldw);
            b1map(m,n,p,2) = 0.15 + exp(-i*pi/10) * exp(-((m-35).^2 + (n-5).^2 + (p-35).^2 )/ fieldw);
            
            b1map(m,n,p,3) = 0.15 + exp(-i*pi/11) * exp(-((m-60).^2 + (n-35).^2 + (p-20).^2 )/ fieldw);
            b1map(m,n,p,4) = 0.15 + exp(-i*pi/12) * exp(-((m-5).^2 + (n-35).^2 + (p-20).^2 )/ fieldw);
            
            b1map(m,n,p,5) = 0.15 + exp(-i*pi/11) * exp(-((m-60).^2 + (n-35).^2 + (p-45).^2 )/ fieldw);
            b1map(m,n,p,6) = 0.15 + exp(-i*pi/12) * exp(-((m-5).^2 + (n-35).^2 + (p-45).^2 )/ fieldw);
            
            b1map(m,n,p,7) = 0.15 + exp(-i*pi/7) * exp(-((m-35).^2 + (n-35).^2 + (p-5).^2 )/ fieldw);
            b1map(m,n,p,8) = 0.15 + exp(-i*pi/8) * exp(-((m-35).^2 + (n-35).^2 + (p-60).^2 )/ fieldw);
        end
    end
end

% Next - make the sensitivities significantly different from each other
% for troubleshooting purposes.
b1map = b1map *10;
%
b1map(:,:,:,2) = 2 * b1map(:,:,:,2);
b1map(:,:,:,4) = 2 * b1map(:,:,:,4);
b1map(:,:,:,6) = 2 * b1map(:,:,:,6);
%}

% make a synthetic object
im = phantom3d(MYDIM);
%im(52:60, 42:60, 42:60) = 2;

% using noise to give more texture to the image
msk = ones(size(im));
msk(im==0) = 0;
mynoise = (smooth3(randn(size(im))));
im = abs(msk .* im.*(2+mynoise));
lightbox(im); colormap parula
drawnow
im = im(:);

% weigh the object by individual coil sensitivity maps 
im_sens = zeros(length(im), Ncoils);
figure
for c = 1:Ncoils
    b1 = b1map(:,:,:,c);  b1 = b1(:);
    im_sens(:,c) = im.*b1;

    subplot(3,3,c)
    b=(reshape(abs(im_sens(:,c)), MYDIM, MYDIM, MYDIM) ); 
    lightbox((b(: ,:, 1:3:end)));
    %caxis([0 2])
    title(sprintf('coil %d', c))
end
colormap parula
drawnow

%%
% image space grid
FOV = 20;
[rx , ry, rz ]= meshgrid( ...
    linspace(-FOV/2,FOV/2,MYDIM), ...
    linspace(-FOV/2,FOV/2,MYDIM), ...
    linspace(-FOV/2,FOV/2,MYDIM));
rx = rx(:);
ry = ry(:);
rz = rz(:);
rlocs = [rx ry rz];

% compute the signal equation for each coil 
Nvox = MYDIM^3;
Nt = length(ks);
signal = zeros(Nt, Ncoils);


fprintf('\nSynthesizing signal for each coil from the signal equation .. ');
fprintf('\nBreaking trajectory into segments for efficiency ...');

tic  

% figure out how many segments to break signal into:
klen = length(kxs);
Nsegments = 40;
while mod(klen, Nsegments)
    Nsegments = Nsegments+1;
end
% size of each segment
kseglen = klen/Nsegments;
inds = zeros(kseglen,1);
tmp = zeros(kseglen, Ncoils);
    

for n=1:Nsegments
    t = tic;
    inds = [1:kseglen] + (n-1)*kseglen;
    ktraj = [kxs(inds), kys(inds), kzs(inds)];
    
    signal(inds, :) = mr_signal(im_sens , ktraj,  rlocs);
    
    fprintf('\n segment %d of %d ... %f seconds', n, Nsegments, toc(t));
end

figure
plot(abs(signal(1:1e4,:)))  

return

%%
function signal = mr_signal(rho, ks,  r)
% function signal = mr_signal(rhomap, ks,  r)
%
% computes the signal equation integral 
% Sum signals over all voxels
%
%   rho : spin density map, weighted by coils  (Nlocs x Ncoils)
%   ks :  kspace trajetory segment (kseglen x 3)
%   r :   image space coordinates   (Nlocs x 3)
%
% calculate the phase imparted by the k space trajecotry at each position
% Then do the weighted sum of the signals from all voxels (weighted by the
% density)

Ncoils = size(rho,2);
kseglen = size(ks,1);
signal = zeros(kseglen, Ncoils);

kphase = exp(-i *2*pi* (ks*r')) ;
for c=1:Ncoils
    signal(:,c)=  kphase*rho(:,c);
end

return