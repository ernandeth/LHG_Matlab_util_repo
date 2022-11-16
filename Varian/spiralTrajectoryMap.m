%   Program which can map out and compare the angle between the designed spiral
%   trajectory and what is measured.

clear; close all;

% Set Constants
gamma = 4257;  %radians/s/G
dt = 4e-6;     % Gradient sample time is seconds
close all;
td = '/home/histo/vnmrsys/data/s_2013020901/';
spDir = 'spiral_01.fid';
KxDir = 'KxMap_01.fid';
KyDir = 'KyMap_01.fid';


% Get Spiral Info
[dta, procpar] = getFid([td '/' spDir]);
% [dta, fheader, bheader] = readfid('fid', -1);
sw = procpar.sw;
ts = (1:procpar.np/2)/sw;

% get Kspace locations
KyMeas = getTraj([td KyDir]);
KxMeas = getTraj([td KxDir]);
[Gyi, Kyi, tg] = getIdealGrads([td '/' spDir '/Waveforms'], 'Gy.GRD', dt, sw, 1:size(dta,1));
[Gxi, Kxi, tg] = getIdealGrads([td '/' spDir '/Waveforms'], 'Gx.GRD', dt, sw, 1:size(dta,1));

% Find absolute values of each trajectory
KM = KxMeas+1i*KyMeas;
KI = Kxi+1i*Kyi;
KMval = abs(KM);
KIval = abs(KI);
angles = (unwrap(angle(KM))-unwrap(angle(KI)))*180/pi;
plot(ts, angles); title('Phase difference between actual and measured trajectory'); ylabel('Degrees')
