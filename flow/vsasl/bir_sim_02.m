% using Meakin 2013 formulation for BIR-8 pulses
% units are seconds, cm, radians, Gauss
% this version is trying to simplify to BIR 4
clear all

Npoints = 500;
Tseg = 2e-3;  % seconds
t=linspace(0,Tseg,Npoints);
zeta = 15;
kappa = atan(60); % 1.5541;
tankappa = 60;
wmax = 0.389e3 *2*pi;  % rad/sec  
Gmax = 0;   % G/cm
RFmax = 20e-6;  % Tesla
RFmax = RFmax * 1e4 ; % to Gauss

gambar = 42.576e6;     % gamma in Hz/T

vel = 0;  % cm/sec.
T1 = 1.6;  % sec
T2 = 0.150;  % sec

amp1 = tanh(zeta * (1 - (t./Tseg))); 
amp2 = tanh(zeta * (t./Tseg)); 

ampsweep = [amp1 amp2];
ampsweep = ampsweep(2:end);

phase1 = -wmax*Tseg*log( abs(cos(kappa*t/Tseg)) / (kappa*tankappa) );
phase2 = -wmax*Tseg*log( abs(cos(kappa*(t/Tseg -1))) / (kappa*tankappa) );

dt = Tseg/Npoints;
frqsweep = diff([phase2 phase1])/dt/(2*pi);

%phase1 = mod(phase1, 2*pi)-pi;
%phase2 = mod(phase2, 2*pi)-pi;

rf1 = RFmax*amp1.*exp(i*phase1);
rf2 = RFmax*amp2.*exp(i*phase2);


figure(1)
subplot(321), plot(amp1)
subplot(322), plot(amp2)
subplot(323), plot(phase1)
subplot(324), plot(phase2)
subplot(313), plot(frqsweep)

figure(2)

plot(ampsweep, frqsweep);

% trap1 = ones(1, Npoints);
% trap1(1:Npoints/4) = linspace(0,1,Npoints/4);
% trap1(Npoints*3/4+1:end) = linspace(1,0,Npoints/4);
% trap1 = trap1 *Gmax;

trap1 = ones(1, 2*Npoints);
trap1(1:Npoints/4) = linspace(0,1, Npoints/4);
trap1(Npoints + Npoints*3/4+1:end) = linspace(1,0,Npoints/4);
trap1 = trap1 *Gmax;

gap = zeros(1, Npoints);
gap2 = zeros(1, 2*Npoints);

% if you want a different flip angle, you do this:
% rf1b = rf1 *exp(-i*pi/4);
% RF = [rf1b gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];

%RF = [rf1 gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];
%Gz = [gap trap1 gap gap (-trap1) gap gap (-trap1) gap gap trap1 gap];
%t = linspace(0, 12*Tseg, 12*Npoints);
% 
RF = [rf1 gap gap  rf2 rf1 gap gap  rf2 rf1 gap gap rf2 rf1 gap gap rf2];
Gz = [gap trap1 gap gap gap (-trap1) gap gap gap gap (-trap1) gap gap gap trap1 gap];
t = linspace(0, 16*Tseg, 16*Npoints);
% 
% 
% RF = [rf1 gap2    gap2 rf2 rf1 gap2   gap2 rf2];
% Gz = [gap trap1   gap2 gap gap  gap2  trap1 gap];

%RF = [gap gap2    gap2 rf2 rf1 gap2   gap2 gap];
%Gz = [gap trap1   gap2 gap gap  gap2  trap1 gap];

duration = length(Gz)*Tseg/Npoints;
t = linspace(0, duration, length(Gz));

zpos = ones(size(Gz)) .* (0 + vel*t);

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

%% unit conversionsl  Blochsim works in Tesla and mili seconds

T1 = T1 * 1e3;      % s to ms
T2 = T2 * 1e3;      % s to ms
dt = dt * 1e3;      % s to ms
beff = beff*1e-4; % G to Tesla

M = blochsim(Mi, beff, T1, T2, dt, nstep);

M = M(1:length(t),:);

close all

%{
for n=1:length(t), 
    cla
    line([0 M(n,1)], [0 M(n,2)], [0 M(n,3)] ), 
    axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
    %line([0 M(n,2)], [0 M(n,3)] ), 
    %axis([ -1 1 -1 1]);
     %
    grid on
    
   
    drawnow, 
end
%}

subplot(313)
plot(t, M(:,3))
grid on
M(end,:)
Mxy = complex(M(:,1), M(:,2));
Mphase = unwrap(angle(Mxy));
fprintf('\nPhase of Mxy = %f rads if traveling at %f cm/s', Mphase(end), vel);

gamma = 2*pi*4.257e3; % rad/s/gauss
Gz2 = Gz;
Gz2(end/2:end) = - Gz2(end/2:end);
phs = gamma* sum(Gz2.*zpos.*dt);
fprintf('\nAlt. Phase of Mxy = %f rads if traveling at %f cm/s', phs, vel);

return
%% The next chunk makes RF pulses for the GE scanner.

dacmax = hex2dec('7ffe');

fp = fopen('BIR_odd.bin','wb','b');
tmp = abs(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_odd.rho BIR_odd.bin

fp = fopen('BIR_odd.bin','wb','b');
tmp = angle(rf1)  ;
tmp = tmp * dacmax / pi;
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_odd.theta BIR_odd.bin

fp = fopen('BIR_even.bin','wb','b');
tmp = abs(rf2);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_even.rho BIR_even.bin

fp = fopen('BIR_even.bin','wb','b');
tmp = angle(rf2);
tmp = tmp * dacmax / pi;
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_even.theta BIR_even.bin

fp = fopen('BIR_whole.bin','wb','b');
tmp = abs([rf2 rf1]) ;
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
fwrite(fp, tmp, 'int16');
fclose(fp);
!xlatebin -o BIR_whole.rho BIR_whole.bin

fp = fopen('BIR_whole.bin','wb','b');
tmp = angle([rf2 rf1]);
tmp = tmp * dacmax / pi;
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

