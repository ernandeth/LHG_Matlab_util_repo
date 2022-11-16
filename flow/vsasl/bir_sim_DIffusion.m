% using Meakin 2013 formulation for BIR-8 pulses
% units are seconds, cm, radians, Gauss
% this version is trying to simplify to BIR 4
% clear all
% close all

doScannerPulses = 1;
doBatch=0
controlcase = 0;
sampTime=4e-6;
Tseg = 3e-3;  % seconds - % used to be 2 ms
%Tseg = 2e-3;  % LHG 4/14/22
Npoints = round(Tseg/sampTime);
t=linspace(0,Tseg,Npoints);
zeta = 15;

kappa = atan(60); % 1.5541;
tankappa = 60;

% kappa = atan(30); % 1.5541;
% tankappa = 30;

% LHG 4/14/22
%wmax =  0.389e3 *2*pi;  % rad/sec  
wmax =  0.35e3 *2*pi;  % rad/sec  
wmax =  0.30e3 *2*pi;  % rad/sec  

Gmax = 1.5;   % G/cm
RFmax = 15e-6;  % Tesla
RFmax = 10e-6;  % Tesla
RFmax = RFmax * 1e4 ; % to Gauss

% messing with it trying to understand:
% tankappa = tankappa*6;
% kappa = atan(tankappa);
% RFmax = RFmax * 0.8;

zpos0 = 0;
if doBatch==0
    vel = 0;  % cm/sec.
end

T1 = 1.6;  % sec
T2 = 0.150;  % sec

% LOW field :  0.55T 
% from doi: 10.1148/radiol.2019190452
T1 = 1.122;
T2 = 0.263;

amp1 = tanh(zeta * (1 - (t./Tseg))); 
amp2 = tanh(zeta * (t./Tseg)); 

phase1 = -wmax*Tseg*log( abs(cos(kappa*t/Tseg)) / (kappa*tankappa) );
phase2 = -wmax*Tseg*log( abs(cos(kappa*(t/Tseg -1))) / (kappa*tankappa) );


dt = Tseg/Npoints;
frqsweep = diff([phase2 phase1])/dt/(2*pi);

% phase1 = mod(phase1, 2*pi)-pi;
% phase2 = mod(phase2, 2*pi)-pi;

rf1 = RFmax*amp1.*exp(i*phase1);
rf2 = RFmax*amp2.*exp(i*phase2);

% mysech = genSech180(1,4);
% rf1 = RFmax * mysech(1:4:end/2);
% rf2 = RFmax * mysech(end/2+1:4:end);

if doBatch==0
figure(1)
subplot(321), plot(amp1)
subplot(322), plot(amp2)
subplot(323), plot(phase1)
subplot(324), plot(phase2)
subplot(313), plot(frqsweep)
end

% 1 ms  gradients at 8 usec resolution
Gseg = 3e-3;  % seconds
Gseg = 8e-3;  % seconds : LHG 4/14/22 - use for diffusion
ramplen = round(0.25*1e-3/sampTime);
NGpoints = round(Gseg/sampTime);

trap1 = ones(1, NGpoints);
ramp = linspace(0,1, ramplen);
trap1(1:ramplen) = ramp;
trap1(end-ramplen+1:end) = ramp(end:-1:1);

gspace = zeros(1,300) ;  % a little gap between the RF and the VENC to mititgate eddy currrents


rfgap = zeros(1, NGpoints);
ggap = zeros(1, Npoints);

% if you want a different flip angle, you do this:
% rf1b = rf1 *exp(-i*pi/4);
% RF = [rf1b gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];

% BIR-8
RF = [rf1 ...
    rfgap ...
    gspace ...
    -rf2*exp(i*pi/2) -rf1*exp(i*pi/2) ...
    rfgap ...
    gspace ...
    -rf2];
% RF = [rf1 rfgap rf2 rf1 rfgap rf2 rf1 rfgap rf2 rf1 rfgap rf2];


Gz = [ggap trap1     gspace ggap ggap  trap1     gspace ggap];
t = [0:length(Gz)-1] * sampTime;

if controlcase
    Gz=0;
end

% BIR 4 
% RF = [rf1    gap2 rf2 rf1 gap2    rf2];
% Gz = [gap    gap2 gap gap  gap2   gap];

% BIR 4 with gradients
%RF = [gap gap2    gap2 rf2 rf1 gap2   gap2 gap];
%Gz = [gap trap1   gap2 gap gap  gap2  trap1 gap];

%duration = length(Gz)*Tseg/Npoints;
%t = linspace(0, duration, length(Gz));

zpos = ones(size(Gz)) .* (zpos0 + vel*t);

if doBatch==0
%%
figure(2)
subplot(311)
area(t, abs(RF)*1000);
grid off
ylabel('B1 amplitude (mG)')
axis tight
set(gca, 'box','off')

subplot(312)
plot(t, angle(RF),'k'); 
grid off
ylabel('B1 phase (rad)')
axis tight
set(gca, 'box','off')

subplot(313)
area(t, Gz)
grid off
ylabel('Grad.(G/cm)')
axis tight
set(gca, 'box','off')
%%
if doScannerPulses
    genScannerPulses(RF, Gz, sampTime*1e3);
end

end


% Put the whole thing together and run the simulation

Mi = [0 0 1];

Bx = real(RF);
By = imag(RF);
Bz = Gz.*zpos;

% add off-resonance :
% Bz = Bz + 500*2*pi/gamma;

nstep = length(RF);
dt = Tseg/Npoints;

beff = [Bx' By'  Bz'];

%% unit conversionsl  Blochsim works in Tesla and mili seconds

T1 = T1 * 1e3;      % s to ms
T2 = T2 * 1e3;      % s to ms
dt = dt * 1e3;      % s to ms
beff = beff*1e-4; % G to Tesla

M = blochsim(Mi, beff, T1, T2, dt, nstep);

M = M(1:length(RF),:);

%{
figure
for n=1:length(t), 
    cla
%    line([0 M(n,1)], [0 M(n,2)], [0 M(n,3)] ), 
    quiver3(0, 0, 0, M(n,1), M(n,2), M(n,3) ), 
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
    %line([0 M(n,2)], [0 M(n,3)] ), 
    %axis([ -1 1 -1 1]);
     %
    grid on
    
   
    drawnow, 
end
%
figure(4)
subplot(311)
plot(t, M(:,1))
grid on
subplot(312)
plot(t, M(:,2))
grid on
subplot(313)
plot(t, M(:,3))
grid on
%}

M(end,:)
Mxy = complex(M(:,1), M(:,2));
Mphase = (angle(Mxy));
fprintf('\nPhase of Mxy = %f rads if traveling at %f cm/s', Mphase(end), vel);

gamma = 2*pi*4.257e3; % rad/s/gauss
Gz2 = Gz;
Gz2(end/2:end) = - Gz2(end/2:end);
phs = gamma* sum(Gz2.*zpos.*dt);
fprintf('\nAlt. Phase of Mxy = %f rads if traveling at %f cm/s', phs, vel);



return
%
%% The next chunk makes RF pulses for the GE scanner.
%
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
%%
% 32 ms. pulse sampled at 8 us
NPOINTS = length(RF);
dacmax = hex2dec('7ffe');

tmp = abs(RF)';
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
filename = sprintf('myVSI_%d.rho.txt', NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

tmp = angle(RF)';
tmp = tmp * dacmax / pi;
tmp = 2*round(tmp/2);
filename = sprintf('myVSI_%d.theta.txt', NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

tmp = Gz';
tmp = tmp * dacmax / max(abs(tmp));
tmp = 2*round(tmp/2);
filename = sprintf('myVSI_%d.grad.txt', NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

%%
% Make a 180 for checking
RF = [rf2 rf1 rfgap];
Gz = [ggap ggap trap1];
NPOINTS = length(RF);
dacmax = hex2dec('7ffe');

tmp = abs(RF)';
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
filename = sprintf('BIR1_train.rho.txt');
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

tmp = angle(RF)';
tmp = tmp * dacmax / pi;
tmp = 2*round(tmp/2);
filename = sprintf('BIR1_train.theta.txt');
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

tmp = Gz';
tmp = tmp * dacmax / max(abs(tmp));
tmp = 2*round(tmp/2);
%tmp(:) =0;
filename = sprintf('BIR1_train.grad.txt');
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);
