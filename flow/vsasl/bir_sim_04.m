% using Meakin 2013 formulation for BIR-8 pulses
% units are seconds, cm, radians, Gauss
% this version is trying to use the BIR 180 pulse in order to get adiabatic
% inversion pulses that are velocity selective

%clear all
%close all
showMag = 1;

Npoints = 500;  % NUmber of points in original 90 pulse
Tseg = 4e-3;  % seconds  : duration of the original 90 pulse
t = linspace(0,Tseg,Npoints);
zeta = 15;
kappa = atan(60); % 1.5541;
tankappa = 60;
wmax = 0.389e3 *2*pi;  % rad/sec  

Gmax = 0.5;   % G/cm
RFmax = 50e-6;  % Tesla
RFmax = RFmax * 1e4 ; % to Gauss

pos0 = 0 ; % cm
%vel = 40 ; % cm/sec
accel = 0;

vel_target = 40;

T1 = 1.6;  % sec
T2 = 0.180;  % sec

% The usual constants
gambar = 42576;     % gamma in kHz/T
gam = gambar*2*pi;  % radians / (ms *Tesla)
mygam = gam / 1000; % radians / (s*Tesla)

amp1 = tanh(zeta * (1 - (t./Tseg))); 
amp2 = tanh(zeta * (t./Tseg)); 

phase1 = -wmax*Tseg*log( abs(cos(kappa*t/Tseg)) / (kappa*tankappa) );
phase2 = -wmax*Tseg*log( abs(cos(kappa*(t/Tseg -1))) / (kappa*tankappa) );


dt = Tseg/Npoints;
frqsweep = diff([phase2 phase1])/dt/(2*pi);

%phase1 = mod(phase1, 2*pi)-pi;
%phase2 = mod(phase2, 2*pi)-pi;

rf1 = RFmax*amp1.*exp(i*phase1);
rf2 = RFmax*amp2.*exp(i*phase2);

rf180 = [rf2 rf1];

figure(1)

subplot(321), plot(amp2)
subplot(322), plot(amp1)
subplot(323), plot(phase2)
subplot(324), plot(phase1)
subplot(313), plot(frqsweep)

trapLength = Npoints;
trap = ones(1, trapLength);
trap(1:trapLength/2) = linspace(0,1, trapLength/2);
trap(trapLength/2+1:end) = linspace(1,0,trapLength/2);
trap = trap *Gmax;

Nsegs=4;
Npoints = length(rf180);
B1segLength = Npoints/Nsegs;
B1module = zeros(1, B1segLength);
Gzmodule = [B1module, trap , -trap];

% This is the phase of the spins traveling at the target velocity
% during one of the Gzmodules.
% We will use this phase as compensation during RF time.
% first calculate the position if it moves with the target velocity and acceleration:
ttmp = [0:length(Gzmodule)-1]*dt;  % a time vector in seconds.
zpos_target = vel_target * ttmp ;

% Then the phase acquired during the module.
Mxy_phase = -gam * cumsum( zpos_target .* Gzmodule)*dt   ;
RF_phase = Mxy_phase(end) ;

B1module = zeros(size(Gzmodule));
B1pad = zeros(1, length(trap)*2 );

gap = zeros(1, Npoints);
gap2 = zeros(1, 2*Npoints);

beg = 1;
fin = round(Npoints/Nsegs);

for n = 1: Nsegs
    
    B1seg = rf180(beg:fin);
    
    % Each pulse will have some extra phase to compensate for the phase
    % acquisred by moving spins in the presence of a gradient.
    B1seg = B1seg * exp(i*RF_phase*(n-1)); 
    
    beg = fin + 1;
    fin = beg + B1segLength-1;
    
    % update the excitation module
    B1module = [B1seg , B1pad];

    modlength = length(B1module);
    
    B1( (n-1)*modlength + 1 : n*modlength) = B1module;
    Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule;
    

end

% B1 = rf180;
% Gz = zeros(size(B1));

% if you want a different flip angle, you do this:
% rf1b = rf1 *exp(-i*pi/4);
% RF = [rf1b gap rf2 rf1 gap rf2 rf1 gap rf2 rf1 gap rf2];

t = linspace(0, 5*Tseg,length(Gz));


% Now the movement of the particle
% calculate the actual position of particle in one module
zpos = pos0 + vel*t + accel*(t.^2);

% position of particle in cm over the whole period
dx = 0;
dy = 0;
dz = zpos;

   
Gx = zeros(size(Gz));
Gy = zeros(size(Gz));

% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz = (dx.*Gx + dy.*Gy + dz.*Gz) ;

    
% off resonance!!
% Bz = Bz + randn(size(Bz))*2e-5;

% these are the transverse fields  due to B1 pulse (in Tesla )
Bx = real(B1);
By = imag(B1);

%%

% Put the whole thing together and run the simulation
beff = [Bx' By'  Bz'];
NSTEPS = length(Bx);

Mi = [0,0,1];
M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
%M = M(1:end-1,:);

%%
% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)



% Make nice plots here
if showMag
    figure(2)
    subplot(211)
    area(t, Gz);
    subplot(212)
    plot(t, Bx,'r'); hold on;
    plot(t, By,'g'); 
    plot(t, abs(B1), 'k'); hold off;
    
    figure(3)
    set(gcf,'Name', 'Magnetization (x,y,z)')
    set(gcf,'Position',[ 600 100 500 600 ])
    timevec = [0:(length(M)-1)]*dt;
    subplot(411) , plot(t, M(:,1))
    subplot(412) , plot(t, M(:,2))
    subplot(413) , plot(t, M(:,3))
    axis([0 t(end) -1 1]);
    subplot(414), plot(t,zpos);
    plot(t, mygam * cumsum(Bz * dt));
    
    grid on
    %
    %     figure(4)
    %     plot(t,dz)
    %     title('Position of particle in time')
    %
    
end


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

