function [signal b1map ] = sim_kspace3d_sense(kxs, kys, kzs, densmap)
% function [signal b1maps ] = sim_kspace3d_sense(kxs, kys, kzs, densmap)
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
MYDIM = 64;
if nargin==4
    MYDIM = size(densmap,1)
end

Ncoils = 8;
ks = [kxs, kys, kzs];
Kmax = max(abs(ks(:)));

% Make some very simplistic complex sensitivity maps
b1map = zeros(MYDIM, MYDIM, MYDIM, Ncoils);
kim = b1map;
fieldw = 5000;
offset = MYDIM/3;

dxc = [1 -1  1  1  1  1 -1 -1];
dyc=  [1  1 -1  1  1 -1 -1 -1];
dzc = [1  1  1 -1 -1 -1  1 -1];

for c=1:8

    dx = MYDIM/2 - offset*dxc(c);
    dy = MYDIM/2 - offset*dyc(c);
    dz = MYDIM/2 - offset*dzc(c);
    
    for m=1:MYDIM
        for n=1:MYDIM
            for p=1:MYDIM
                b1map(m,n,p,c) = exp(-i*c*pi/9 ) * ...
                    exp(-((m-dx).^2 + (n-dy).^2 + (p-dz).^2 )/ fieldw);
            end
        end
    end
end
% Next - make the sensitivities significantly different from each other
% for troubleshooting purposes.
%{
b1map(:,:,:,2) = 2 * b1map(:,:,:,2);
b1map(:,:,:,4) = 2 * b1map(:,:,:,4);
b1map(:,:,:,6) = 2 * b1map(:,:,:,6);
%}

% Normalize the sense maps here abd (Ssquaress) should be 1
btmp = sum(abs(b1map).^2, 4);
orthoview(abs(b1map(:,:,:,4)))

if nargin==4
    im = densmap;
else
    % make a synthetic object
    im = phantom3d(MYDIM);
    %im(52:60, 42:60, 42:60) = 2;

    % using noise to give more texture to the image
    msk = ones(size(im));
    msk(im==0) = 0;
    mynoise = (smooth3(randn(size(im))));
    im = abs(msk .* im.*(2+mynoise));
end

lightbox(im); colormap parula
drawnow
im = im(:);

% weigh the object by individual coil sensitivity maps 
im_sens = zeros(length(im), Ncoils);
figure
b_rms = 0;
for c = 1:Ncoils
    b1 = b1map(:,:,:,c);  
    b1 = b1(:);
    im_sens(:,c) = im.*b1;

    subplot(3,3,c)
    b=(reshape((im_sens(:,c)), MYDIM, MYDIM, MYDIM) ); 
    orthoview(abs(b));
    %caxis([0 2])
    title(sprintf('coil %d', c))

    b_rms = b_rms + b.^2;
end
b_rms = sqrt(b_rms);
colormap parula
drawnow

figure
lbview(abs(b_rms))
title('RMS sum of coils')
%%
% image space grid
FOV = 24;
[rx , ry, rz ]= meshgrid( ...
    linspace(-FOV/2,FOV/2,MYDIM), ...
    linspace(-FOV/2,FOV/2,MYDIM), ...
    linspace(-FOV/2,FOV/2,MYDIM));

rx = rx(:);
ry = ry(:);
rz = rz(:);


% compute the signal equation for each coil 
Nvox = MYDIM^3;
Nt = length(ks);
signal = dlarray(zeros(Nt, Ncoils));


fprintf('\nSynthesizing signal for each coil (signal equation) ... ');
%fprintf('\nBreaking trajectory into segments for efficiency ...');
 
gpuDevice(1)
ktraj = gpuArray([kxs, kys,  kzs]);
rlocs = gpuArray([rx ry rz]);
im_sens = gpuArray(im_sens);

for n=1:Ncoils
    
    signal(:, n) = mr_signal(im_sens(:,n) , ktraj,  rlocs);
    
    
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
%   rho : spin density map, weighted by the coil  (Nlocs x 1)
%   ks :  kspace trajetory segment (kseglen x 3)
%   r :   image space coordinates   (Nlocs x 3)
%
% calculate the phase imparted by the k space trajecotry at each position
% Then do the weighted sum of the signals from all voxels (weighted by the
% density)


klen = size(ks,1);
signal = zeros(klen, 1);

%----
% figure out how many segments to break signal into:
% each segment will have a fixed size
kseglen = 100;
Nsegments = floor(klen/kseglen)+1;
kseglen_end = mod(klen , kseglen);  % leftover points for last segment
inds = zeros(kseglen,1);
tmp = zeros(kseglen, 1);

    %kphase = exp(-i *2*pi* (ks*r')) ;
%
for s=1:Nsegments-1
    %fprintf('\rSegment %d of %d. ', s, Nsegments);
    %tic
    inds = [1:kseglen] + (s-1)*kseglen;
    kphase = exp(-i *2*pi* (ks(inds,:)*r')) ;
    signal(inds)=  kphase*rho;
    %toc 
end

% now the remainder of the points ...
inds = [1:kseglen_end] + (Nsegments-1)*kseglen;
kphase = exp(-i *2*pi* (ks(inds,:)*r')) ;
%}
    signal(inds)=  kphase*rho;
%---------


%{
kphase = exp(-i *2*pi* (ks*r')) ;
for c=1:Ncoils
    signal(:,c)=  kphase*rho(:,c);
end
%}
return