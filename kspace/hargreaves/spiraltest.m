% testing spirals
gambar = 4.257e1;     % kHz/mT

% uniform spiral grad
opfov = 20; % cm
fov = opfov/100; % in m
opxres = 64;  % matrix size
gts = 4e-3;  % gradient sampling time (ms)
gslew = 150; % mT/m/ms
gamp = 40; % mT/m
nl = 1;  % interleaves
% conversions for strange units of spiralgradlx3
[g,k,t,s]=spiralgradlx3(fov*100,opxres,gts/1000,gslew,gamp/10,nl);
% convert g to mT/m
g = g*10;
% ramp down
lramp = ceil(abs(g(end))/gslew/gts);
gnew = [g ((lramp-1):-1:0)/lramp*g(end)];
% rephase
garea = sum(gnew)*gts;
gnew = [gnew dotrap(-garea,gslew,gamp,gts)];

% do spiral in/out
halfl = length(gnew);
gnew = [(gnew(end:-1:1)) gnew];

% calculate k-space and slew
time = [1:length(gnew)]*gts;
knew = cumsum(gnew)*gts*gambar; %1/cm
knew = knew*fov; % samples of 1/FOV
snew = [0 diff(gnew)]/gts;

% plotting
subplot([311])
plot(time,real(gnew),time,imag(gnew))
title('grads')
subplot([312])
plot(time,real(snew),time,imag(snew))
title('slew')
subplot([337])
plot(knew)
subplot([337])
plot(knew(1:halfl))
hold on; plot(knew(halfl+1:end),'r'); hold off
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')
title('spiral in-out')
subplot([338])
plot(knew(1:halfl))
title('spiral in')
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')
subplot([339])
plot(knew(halfl+1:end),'r')
title('spiral out')
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')

%%
% now try Hargreaves VD spiral version
% code, including c, mex files, etc. is here
% https://mrsrl.stanford.edu/~brian/vdspiral/
%

vdfactor = 1/3;
Fcoeff = [opfov -opfov*(1-vdfactor)] 	% FOV decreases linearly from fov to fov*vdfactor
rmax = (opxres/opfov)/2;		% max kr cm^(-1)

disp('Calculating Gradient');
figure
% again convert to strange units
[k,g,s,time,r,theta] = vds(gslew*100,gamp/10,gts/1000,nl,Fcoeff,rmax);
g = g*10;

% ramp down
lramp = ceil(abs(g(end))/gslew/gts);
gnew = [g ((lramp-1):-1:0)/lramp*g(end)];
% rephase
garea = sum(gnew)*gts;
gnew = [gnew dotrap(-garea,gslew,gamp,gts)];

% do spiral in/out
halfl = length(gnew);
gnew = [(gnew(end:-1:1)) gnew];

% calculate k-space and slew
time = [1:length(gnew)]*gts;
knew = cumsum(gnew)*gts*gambar; %1/cm
knew = knew*fov; % samples of 1/FOV
snew = [0 diff(gnew)]/gts;

% plotting
subplot([311])
plot(time,real(gnew),time,imag(gnew))
title('grads')
subplot([312])
plot(time,real(snew),time,imag(snew))
title('slew')
subplot([337])
plot(knew)
subplot([337])
plot(knew(1:halfl))
hold on; plot(knew(halfl+1:end),'r'); hold off
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')
title('spiral in-out')
subplot([338])
plot(knew(1:halfl))
title('spiral in')
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')
subplot([339])
plot(knew(halfl+1:end),'r')
title('spiral out')
axis( [-opxres/2 opxres/2 -opxres/2 opxres/2]*1.2); axis('square')