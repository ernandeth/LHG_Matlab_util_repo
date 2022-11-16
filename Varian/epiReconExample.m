%  Matlab script as example for recon algorithms.


% setPathRecon;   %  Set Path
% setup;          % Setup Image Reconstruction Toolbox
clear;
close all;

topdir = '/home/histo/vnmrsys/data/s_2013040901';
spscan = 'spiral_01.fid';
KxTraj = 'KxMap_02.fid';
KyTraj = 'KyMap_02.fid';


% Matlab based recon.
im = reconstructSpiralEpi(topdir,spscan, 20,KxTraj,KyTraj);
figure(1),  imagesc(fliplr(abs(im(32:96,32:96,15)))); colormap(gray), colorbar, title('Matlab Regrid')

% Fesslerian Nufft recon
% picture = nufftRecon(topdir,spscan,KxTraj,KyTraj);
%  figure(2), imagesc(abs(imrotate(picture,-90))); colormap(gray); colorbar; title('Fesslerian Reconstruction');