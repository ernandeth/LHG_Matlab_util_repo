%  Matlab script as example for recon algorithms.


% setPathRecon;   %  Set Path
% setup;          % Setup Image Reconstruction Toolbox
clear;
% close all;

topdir = '/home/histo/vnmrsys/data/s_2014090501';
spscan = 'spiral_02.fid';
KPeTraj = 'KPeMap_02.fid';
KRoTraj = 'KRoMap_02.fid';
text = 'Spiral 1 Shots 205 in 205';
filename = 'spiral_1Shot_305in205.jpg';

% Matlab based recon.
 im = reconstructSpirals(topdir,spscan, KRoTraj,KPeTraj);
[a b] = size(im);
traj = figure(1);  imagesc((abs(im(a/4:end-a/4, b/4:end-b/4)))); colormap(gray), colorbar, title([text ': With Mapping'])

%im = reconstructSpirals(topdir,spscan);%, KxTraj,KyTraj);
%[a b] = size(im);
%noTraj = figure(2);  imagesc(fliplr(abs(im(a/4:end-a/4, b/4:end-b/4)))); colormap(gray), colorbar, title([text ': No Mapping'])
display('Recon Done');

% Print file to jpeg
%print(traj, '-djpeg', filename);
% % Fesslerian Nufft recon
%  picture = nufftRecon(topdir,spscan,KxTraj,KyTraj);
%   figure(3), imagesc(abs(imrotate(picture,-90))); colormap(gray); colorbar; title('Fesslerian Reconstruction');
