% using Meakin 2013 formulation for BIR-8 pulses
% units are seconds, cm, radians, Gauss
% this version is trying to simplify to BIR 4
clear all
close all
clc

Npoints = 250;
Tseg = 2e-3;  % seconds
Tgap = 100e-6;
t=linspace(0,Tseg,Npoints);

gamma = 2*pi*4.257e3; % rad/s/gauss

zeta = 15;
kappa = atan(60); % 1.5541;
tankappa = 60;
wmax = 0.389e3 *2*pi;  % rad/sec  
Gmax = 2;   % G/cm
RFmax = 20e-6;  % Tesla
RFmax = RFmax * 1e4 ; % to Gauss

venc = 2; % cm/sec.
Gmax = pi/(2*gamma * venc * (2*2*Tseg+2*Tseg)*(Tseg*3/4));
vel = 2;  % cm/sec.
T1 = 1.6;  % sec
T2 = 0.150;  % sec

amp1 = tanh(zeta * (1 - (t./Tseg))); 
amp2 = tanh(zeta * (t./Tseg)); 

phase1 = -wmax*Tseg*log( abs(cos(kappa*t/Tseg)) / (kappa*tankappa) );
phase1 = phase1 - phase1(1);
phase2 = -wmax*Tseg*log( abs(cos(kappa*(t/Tseg -1))) / (kappa*tankappa) );
phase2 = phase2 - phase2(end);

dt = Tseg/Npoints;
frqsweep = diff([phase2 phase1])/dt/(2*pi);

rf1 = RFmax*amp1.*exp(1i*phase1);
rf2 = RFmax*amp2.*exp(1i*phase2);

trap1 = ones(1, Npoints);
trap1(1:Npoints/4) = linspace(0,1, Npoints/4);
trap1(Npoints*3/4+1:end) = linspace(1,0,Npoints/4);
trap1 = trap1 *Gmax;

gap = zeros(1, Npoints);
greatgap = zeros(1, 1*Npoints);
Ngap = floor(Npoints/Tseg*Tgap);
tgap = zeros(1, Ngap);

RF = [gap rf1 gap   rf2 rf1  gap     rf2 rf1 gap      rf2 rf1  gap   rf2 greatgap];
Gz = [gap gap trap1 gap gap (-trap1) gap gap (-trap1) gap gap  trap1 gap greatgap];
t = linspace(0, (13+1)*Tseg, (13+1)*Npoints);

% RF4 = [gap rf1 tgap rf2 rf1 tgap rf2 gap];
% Gz = [gap gap tgap gap gap tgap gap gap];
% t = linspace(0, 6*Tseg+2*Tgap, 6*Npoints+2*Ngap);

zpos = ones(size(Gz)) .* (0 + vel*t);

figure(5)
subplot(311)
plot(t, Gz)
subplot(312)
plot(t, angle(RF),'g'); hold on;
plot(t, abs(RF)*10); hold off;

RF0 = RF;

% Put the whole thing together and run the simulation
n=1;myMz = [];
myRF = 0:0.05:3;

T1 = T1 * 1e3;      % s to ms
T2 = T2 * 1e3;      % s to ms

%%
for ii = myRF
    
    disp('-----------------------------------');
    disp(ii);
    
    RF = RF0 * ii;
    
    Mi = [0 0 1];

    Bx = real(RF);
    By = imag(RF);
    Bz = Gz.*zpos;

    nstep = length(RF);
    dt = Tseg/Npoints;
    dt = dt * 1e3;      % s to ms

    beff = [Bx' By' Bz'];

    %% unit conversionsl  Blochsim works in Tesla and mili seconds

    beff = beff*1e-4; % G to Tesla

    %M = blochsim(Mi, beff, T1, T2, dt, nstep);
    M = myBlochSim(Mi, beff, T1, T2, dt);
    
    M = M(1:length(t),:);

    figure(67);

    subplot(311)
    plot(t, M(:,1))
    grid on
    subplot(312)
    plot(t, M(:,2))
    grid on
    subplot(313)
    plot(t, M(:,3))
    grid on
    drawnow;
    
    pause(0.2),
    myMz(n) = M(end,3);
    n=n+1;
end
figure;
plot(myRF, myMz)

M(end,:)
Mxy = complex(M(:,1), M(:,2));
Mphase = unwrap(angle(Mxy));
fprintf('Phase of Mxy = %f rads if traveling at %f cm/s\n', Mphase(end), vel);

Gz2 = Gz;
Gz2(end/2:end) = - Gz2(end/2:end);
phs = gamma* sum(Gz2.*zpos.*dt);
fprintf('Alt. Phase of Mxy = %f rads if traveling at %f cm/s\n', phs, vel);

return
