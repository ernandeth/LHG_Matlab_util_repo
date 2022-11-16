% using Meakin 2013 formulation for BIR-8 pulses
% units are seconds, radians, Gauss


Npoints = 1000;
Tseg = 1e-3;
t=linspace(0,Tseg,Npoints);
zeta = 15;
kappa = 1.5541;
wmax = 5e3;
Gmax = 1.5;
RFmax = 50e-3;
v = 0;
T1 = 1.6;
T2 = 0.100;

amp1 = tanh(zeta * (1 - (t./Tseg))); 
amp2 = tanh(zeta * (t./Tseg)); 

phase1 = -wmax*Tseg*log( abs(cos(kappa*t/Tseg)) / (kappa*tan(kappa)) );
phase2 = -wmax*Tseg*log( abs(cos(kappa*(t/Tseg -1))) / (kappa*tan(kappa)) );

% phase1 = mod(phase1, 2*pi)-pi;
% phase2 = mod(phase2, 2*pi)-pi;

rf1 = RFmax*amp1.*exp(-i*phase1);
rf2 = RFmax*amp2.*exp(-i*phase2);


figure(1)
subplot(221), plot(amp1)
subplot(222), plot(amp2)
subplot(223), plot(phase1)
subplot(224), plot(phase2)

trap1 = ones(1, Npoints);
trap1(1:Npoints/4) = linspace(0,1,Npoints/4);
trap1(Npoints*3/4+1:end) = linspace(1,0,Npoints/4);
trap1 = trap1 *Gmax;

gap = zeros(1, Npoints);

% if you want a different flip angle, you do this:
% rf1b = rf1 *exp(-i*pi/4);
% RF = [rf1b gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];

RF = [rf1 gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];
Gz = [gap trap1 gap gap (-trap1) gap gap (-trap1) gap gap trap1 gap];
t = linspace(0, 12*Tseg, 12*Npoints);

RF = [rf1 gap gap  rf2 rf1 gap gap  rf2 rf1 gap gap rf2 rf1 gap gap rf2];
Gz = [gap trap1 gap gap gap (-trap1) gap gap gap gap (-trap1) gap gap gap trap1 gap];
t = linspace(0, 16*Tseg, 16*Npoints);

zpos = ones(size(Gz)) .* (5 + v*t);

figure(2)
subplot(311)
plot(t, Gz)
subplot(312)
plot(t, angle(RF),'g'); hold on;
plot(t, abs(RF)*10); hold off;



% Put the whole thing together and run the simulation

Mi = [0 0 1];

Bx = real(RF);
By = imag(RF);
Bz = Gz.*zpos;

nstep = length(RF);
dt = Tseg/Npoints;

beff = [Bx' By'  Bz'];
M = blochsim(Mi, beff, T1, T2, dt, nstep);
%M = M(1:end-1,:);

subplot(313)
plot(t, M(:,3))
M(end,:)

%% The next chunk makes RF pulses for the GE scanner.

wmax = hex2dec('7ffe');

fp = fopen('BIR_odd.bin','wb','b');
tmp = abs(rf1);
tmp = tmp * wmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_odd.rho BIR_odd.bin

fp = fopen('BIR_odd.bin','wb','b');
tmp = angle(rf1)  ;
tmp = tmp * wmax / pi;
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_odd.theta BIR_odd.bin

fp = fopen('BIR_even.bin','wb','b');
tmp = abs(rf2);
tmp = tmp * wmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_even.rho BIR_even.bin

fp = fopen('BIR_even.bin','wb','b');
tmp = angle(rf2);
tmp = tmp * wmax / pi;
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_even.theta BIR_even.bin

fp = fopen('BIR_whole.bin','wb','b');
tmp = abs([rf2 rf1]) ;
tmp = tmp * wmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_whole.rho BIR_whole.bin

fp = fopen('BIR_whole.bin','wb','b');
tmp = angle([rf2 rf1]);
tmp = tmp * wmax / pi;
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_whole.theta BIR_whole.bin


%%
figure(3)
fp=fopen('BIR_even.rho','rb','b');
x=fread(fp,1068,'int16');
fclose(fp)
subplot(231)
plot(x)

subplot(234)
fp=fopen('BIR_even.theta','rb','b');
x=fread(fp,1068,'int16');
fclose(fp)
plot(x)

fp=fopen('BIR_odd.rho','rb','b');
x=fread(fp,1068,'int16');
fclose(fp)
subplot(232)
plot(x)

subplot(235)
fp=fopen('BIR_odd.theta','rb','b');
x=fread(fp,1068,'int16');
fclose(fp)
plot(x)

fp=fopen('BIR_whole.rho','rb','b');
x=fread(fp,2*1068,'int16');
fclose(fp)
subplot(233)
plot(x)

subplot(236)
fp=fopen('BIR_whole.theta','rb','b');
x=fread(fp,2*1068,'int16');
fclose(fp)
plot(x)


% 
% subplot(212)
% fp=fopen('sech_7360.theta','rb','b');
% x=fread(fp,1068,'int16');
% fclose(fp)
% plot(x)
% 
% subplot(212)
% fp=fopen('sech_7360.rho','rb','b');
% x=fread(fp,1068,'int16');
% fclose(fp)
% plot(x)

