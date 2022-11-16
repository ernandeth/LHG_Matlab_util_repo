%  Matlab script that reconstructes spiral scans

clear;
%close all;


controlNum =7;  %number of desired fid file

% Reconstruct the Images
for p = 1:numel(controlNum)
numToScan = controlNum(p);

topdir = '/Users/hernan/data/girf/s_2014111003';
topdir = '/Users/hernan/data/girf/s_2014111003';

spscan = ['spiralHistoGradKernels_' sprintf('%02.0f',numToScan) '.fid'];
KPeTraj = ['KPeHistoMapGradKernels_' sprintf('%02.0f', 1) '.fid'];
KRoTraj = ['KRoHistoMapGradKernels_' sprintf('%02.0f', 1) '.fid'];
text = 'Spiral 1 Shots 205 in 205';
filename = 'spiral_1Shot_305in205.jpg';
girf_file = '/Users/hernan/data/girf/girf_parms_141110.mat';

% [xre xim xmxf] = gen_girf('/Users/hernan/data/girf/s_2014111004/chirpInSpiralCode_03.fid');
% [yre yim ymxf] = gen_girf('/Users/hernan/data/girf/s_2014111004/chirpInSpiralCode_02.fid');
% [zre zim zmxf] = gen_girf('/Users/hernan/data/girf/s_2014111004/chirpInSpiralCode_01.fid');

% save /Users/hernan/data/girf/girf_parms_141110.mat  xre xim xmxf yre yim ymxf zre zim zmxf

% [xre xim xmxf] = gen_girf('/Users/hernan/data/girf/s_2014111001/chirpInSpiralCode_06.fid');
% [yre yim ymxf] = gen_girf('/Users/hernan/data/girf/s_2014111001/chirpInSpiralCode_05.fid');
% [zre zim zmxf] = gen_girf('/Users/hernan/data/girf/s_2014111001/chirpInSpiralCode_04.fid');
% 
% save /Users/hernan/data/girf/girf_parms_141110.mat  xre xim xmxf yre yim ymxf zre zim zmxf


% Matlab based recon.
 [pictureMeasured, pictureCorrected, pictureRaw]  = ...
     reconArraySpiralsWith_Without_GIRF(topdir,spscan, KRoTraj,KPeTraj, girf_file);
 
end
% 
