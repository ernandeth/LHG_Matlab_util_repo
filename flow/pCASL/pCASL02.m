% simulation for pseudoCASL experiment based on blochsim

% this section is skipped to allow for iterations when the doBatch flag is
% turned on.  YOu do this when you are calling the script from a batch
% file, e.g. - see pCASL_batch01.m

doBatch = 0;

if ~doBatch
    flip_ang = 25;
    isTag = 1;
    vel = 30;
    %vel = 0;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    tag_loc = 6.5 ; % cm
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.2;  % this is the gradient imbalance used to gain some phase
    fraction_imbalance = 0.2;
    
    fraction_imbalance = (Gz_ss + Gz_phaser)/(2*Gz_ss);
    
    
    zpos = 0;
    close all
end



gambar = 42570;    % gamma/2pi in kHz/T
gam = gambar*2*pi;

T1 = 1600;  %ms
T2 = 250;   %ms

B0 = 3;    % keep things in Tesla

dt = 1e-3; % 1 us expressed in millliseconds.  dt converts from samples to ms
nstep = 500 / dt;  % 500 ms @  1 us steps.

t = linspace(0,500,nstep)'; % time vector in ms

flip_ang = flip_ang * pi /180;  % flip angle in radians

vel = vel * 1e-6 ; % cm/s scaled to cm/us

t_ramp = t_ramp /dt;

Mi = [0, 0, 1]; % initial condition: tipped into x axis.

% initialize gradients in Tesla/cm.
Gx = zeros(nstep,1);
Gy = zeros(nstep,1);
Gz = zeros(nstep,1);

% initialize the B vector (tesla)
Bx = zeros(nstep,1);
By = zeros(nstep,1);
Bz = zeros(nstep,1);

% Make a quick Gx gradient waveform in G/cm.
% Gx(20/dt:30/dt) = 1;
% Gx(31/dt:40/dt) = -1;

% position in cm
dx = 0;
dy = 0;

dz = -nstep/2:nstep/2-1;
dz = dz' * vel;

if vel==0
    dz = zpos';
end

% stationary case:
%dz(:) = 5;

% calculate transmit frequency offset for desired tag location
% in order to achieve the inversion, 
RF_dw = tag_loc * gambar * Gz_ss * 1e-4; % this is in kHz

% we also need to compute the  gradient moment from the extra gradient area
% in the refocusing pulse (Gz_phaser).  This is what produces the phase accumulation
% the labeling location.  This is a function of the gradient and the
% location of the inversion pulse
% RF_theta = tag_loc * gambar * Gz_phaser * 1e-4 * vel*(0.25)^2 
%RF_theta = tag_loc * gambar * Gz_phaser * 1e-4 * (0.25);
%RF_theta = rem(RF_theta, 2*pi);


% B1 field in mGauss
B1 = zeros(nstep,1);

%%
% %define a 6 milllisecond sinc pulse with two lobes on each side
% rfdt = 4*pi / (6/dt) ;
% mySinc = sinc(-2*pi : rfdt : 2*pi);
% % scale the pulse to 120mGauss max
% mySinc = mySinc * 120;  
% % conversion from mGauss to Tesla
% mySinc = mySinc * 1e-7;
%%

% define a 0.5 ms Hanning pulse and Normalize it so the area is 1.  
myHanning = hanning (0.5 / dt);
myHanning = myHanning / (sum(myHanning).*dt);
% scale it to achieve flip angle - result is in Tesla
% the area required for flip_angle is  = flip / gamma
% and put it in Tesla units
myHanning = myHanning * flip_ang / gam;

BandWidth = 1/0.5; % in KHz for each individual pulse
sl_thick = (BandWidth) / (gambar * Gz_ss );

% shift the frequency of the RF pulse by modulating the RF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRL add equivalent sections to what is in pulse sequence:
my_radians_per_CM = gam * fraction_imbalance*(Gz_ss*1e-4)*(length(myHanning)*dt) ;  %0 + pw_rf1 because no ramp time
% notches occur every 2*PI 
my_notch_distance = 2*pi / (my_radians_per_CM);
% determine which notch is closest to our desired tagging location, then we
% will shift that one by an appropriate amount using a phase increment
% between each RF pulse 
closest_notch = round(tag_loc / my_notch_distance);
RF_theta = (tag_loc - closest_notch*my_notch_distance)/my_notch_distance*2*pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% integrate the pulses into the sequence - in Alsop's abstract: 500us
% pulses every 800 us
cnt=0;

doSign = 1;
for t = t_ramp + 1: RF_spacing/dt : 400/dt
    
    if ~isTag
        doSign = -doSign;
    end

    FMx = cos(RF_dw  * 2*pi * [1:length(myHanning)]' * dt + RF_theta*cnt);
    FMy = sin(RF_dw  * 2*pi * [1:length(myHanning)]' * dt + RF_theta*cnt);
    cnt = cnt+1;

    myHanning_x = myHanning .* FMx;
    myHanning_y = myHanning .* FMy;
    
    % a train of RF pulses
    B1(t + t_ramp : t + t_ramp + length(myHanning) -1 ) = ...
        (myHanning_x  + i * myHanning_y) * doSign;
    
    % slice select gradient
    Gz(t + t_ramp : t + t_ramp + length(myHanning) -1 ) = Gz_ss;
    
    % refocuser for the slice select gradient
    Gz( t + length(myHanning) + 2*t_ramp : ...
        t + length(myHanning)*1.5 + 2*t_ramp ) = ...
        - 2*Gz_ss*(1-fraction_imbalance);
end
    
   
%B1(1100:1100+length(mySinc)-1) = mySinc;

% Throw in an ARTIFACTUAL gradient - B0 inhomogeneity
%Gz = Gz + 0.05;

% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz = (dx.*Gx + dy.*Gy + dz.*Gz) * 1e-4;

% these are the transverse fields - due to B1 pulse.
Bx = real(B1);
By = imag(B1);

% Put the whole thing together and run the simulation
beff = [Bx By  Bz];
M = blochsim(Mi, beff, T1, T2, dt, nstep);
M = M(1:end-1,:);

% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)
phi = gambar * sum(Bz * dt);

if ~doBatch
    
    figure
    subplot(411) , plot(Gx)
    subplot(412) , plot(Gy)
    subplot(413) , plot(Gz)
    subplot(413) , hold on, plot(1e5*real(B1) ,'r')
    subplot(414) , plot(real(B1)), hold on, plot(imag(B1),'r')

    figure

    subplot(211),
    plot(Gz(1:2/dt))
    hold on
    plot(1e5*real(B1(1:2/dt)) ,'r') , plot(1e5*imag(B1(1:2/dt)) ,'g')
    hold off
    dz = dz(1:size(M,1));
    subplot(212), plot(dz,M(:,3))
    %axis([1 length(M) -1.1 1.1])
    grid on
end

%    figure
%    subplot(311) , plot(M(:,1))
%    subplot(312) , plot(M(:,2))
%    subplot(313) , plot(M(:,3))


    
    drawnow


