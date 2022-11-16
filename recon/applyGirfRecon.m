% Script that reads in a K-space trajectory, calculates input gradient
% waveforms, and also some impulse-response function coefficients.  It then
% calculates a impulse respone waveform, and applies it to the input
% gradients.  We then see if the input gradients are anything like the
% measured K-space trajectory.


% setup constants
topDir = '/Users/hernan/data/girf/s_2014082001';
spiralName = 'spiralHistoGradKernels_25.fid';
KPeName = 'KPeMapWithGradKernels_02.fid';
KRoName = 'KRoMapWithGradKernels_02.fid';
coefficNames = 'girf_parms.mat';

% Get coefficients
load([topDir '/' coefficNames]);
Nyq_girf = Nyq;

% Load in Spiral Data
[dtaArray, procpar] = getSpiralFidArray([topDir '/' spiralName]);

%Calculate K-space trajectories
KRoMeas = getTraj([topDir '/' KRoName]);
KPeMeas = getTraj([topDir '/' KPeName]);


% Calculate input functions
[Gpei, Kpei, Groi, Kroi, ts, GPe, KPe, GRo, KRo, tg] = getIdealGrads(procpar);

frequencies = (-numel(ts)/2:(numel(ts)/2-1))/ts(end);


idx = find(abs(frequencies)<=1e4);


Nyq  = max(frequencies); 
Nyq = 1/(2*( ts(2) - ts(1)))

% coronal images:
p = allimagp(3,:);  % z gradient, imaginary component
imjunk = polyval(p,linspace(0,Nyq, length(Gpei)));
% rebuild the negative part of the spectrum
imjunk = [ -imjunk(end:-1:1) imjunk];
    
p = allrealp(3,:);  % z gradient, real component
realjunk = polyval(p,linspace(0,Nyq, length(Gpei)));
% rebuild the negative part of the spectrum
realjunk = [ -realjunk(end:-1:1) realjunk];

% the 'ideal' impulse response function:
Hz_nice = complex(realjunk, imjunk);


Hz_nice(~idx) = 0;
