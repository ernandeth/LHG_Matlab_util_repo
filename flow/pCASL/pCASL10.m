% simulation for pseudoCASL experiment based on blochsim

% this section is skipped to allow for iterations when the doBatch flag is
% turned on.  YOu do this when you are calling the script from a batch
% file, e.g. - see pCASL_batch01.m
%
% this version is intended to test whether we can afford to use smaller
% gradients and still get inversions.

%%
doBatch = 0;

if ~doBatch
    tag_loc = 10;   
    flip_ang = 22.5;
    B1_ave = 2e-6;
    isTag = 1;
    vel = 50;
    slomo =0;
    eta = 0.4;   %
    t_ramp = 0.12;  % ramp time in ms.
    Gz_err = 0;  % artefactual gradient.  accounting for susceptibility
    phase_correction = 0;
    
    % eta = 0.5 means balanced slice-select and refocuser : 
    % the refocuser area is half the ss area
    slomo = 0;
    close all
end
%%

gambar = 42570;    % gamma/2pi in kHz/T
gam = gambar*2*pi;

T1 = 1660;  %ms
T2 = 100;   %ms

B0 = 3;    % keep things in Tesla

dt = 1e-3; % 1 us expressed in millliseconds.  dt converts from samples to ms
nstep = 500 / dt;  % 500 ms @  1 us steps.

% description of slice-select waveforms
RF_spacing = 1 ; % time between RF pulses in ms.
t_ramp = 0.12;  % ramp time in ms.

Gz_ss = 1;  % slice select gradient in G/cm.  <------
t_Gz_ss = 0.5 ; % ms
G_ave = 0.1;
B1_ave = 2e-6; % Tesla

trap_ss = [
    linspace(0,1,t_ramp/dt)';
    ones(t_Gz_ss/dt,1);
    linspace(1, 0,t_ramp/dt)';
    ] ;

trap_ss_r = [
    linspace(0,1,t_ramp/dt)';
    linspace(1, 0,t_ramp/dt)';
    ] ;

Area_ss = Gz_ss*sum(trap_ss);
Gz_ss_r = (G_ave*(length(trap_ss)+length(trap_ss_r)) - Area_ss ) / (t_ramp/dt);
Area_ss_r = -Gz_ss_r*sum(trap_ss_r);

A_xtra = Area_ss_r + Area_ss;


trap_ss = trap_ss * Gz_ss;
trap_ss_r = trap_ss_r * Gz_ss_r;
 

% unit conversions
t = linspace(0,500,nstep)'; % time vector in ms

vel = vel * 1e-6 ; % cm/s scaled to cm/us
t_ramp = t_ramp /dt;
t_Gz_ss = t_Gz_ss/dt;

Mi = [0, 0, 1]; % initial condition: tipped into x axis.

%%
% initialize gradients in Tesla/cm.
Gx = zeros(nstep,1);
Gy = zeros(nstep,1);
Gz = zeros(nstep,1);

% initialize the B vector (tesla)
Bx = zeros(nstep,1);
By = zeros(nstep,1);
Bz = zeros(nstep,1);


% B1 field in mGauss
B1 = zeros(nstep,1);

% Make a quick Gx gradient waveform in G/cm.
% Gx(20/dt:30/dt) = 1;
% Gx(31/dt:40/dt) = -1;

%%
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

%%
% %define a 6 milllisecond sinc pulse with two lobes on each side
% rfdt = 4*pi / (6/dt) ;
% mySinc = sinc(-2*pi : rfdt : 2*pi);
% % scale the pulse to 120mGauss max
% mySinc = mySinc * 120;
% % conversion from mGauss to Tesla
% mySinc = mySinc * 1e-7;
%%

% define a 0.5 ms Hanning pulse and Normalize it so the area is B1_ave.
myHanning = hanning (0.5 / dt);
myHanning = myHanning * B1_ave / sum(myHanning)/dt;
flip_ang = rad2deg(sum(myHanning)*dt*gam)
BandWidth = 1/0.5; % in KHz for each individual pulse
sl_thick = (BandWidth) / (gambar*1e-4 * Gz_ss );


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRL add equivalent sections to what is in pulse sequence:
%{
fraction_imbalance = (Gz_ss + Gz_rephaser + xtra_phase)/(2*Gz_ss);
my_radians_per_CM = (Gz_ss*1e-4)*(length(myHanning)*dt) * gam * fraction_imbalance;
%0 + pw_rf1 because no ramp time
% notches occur every 2*PI
my_notch_distance = 2*pi / (my_radians_per_CM);
% determine which notch is closest to our desired tagging location, then we
% will shift that one by an appropriate amount using a phase increment
% between each RF pulse
closest_notch = round(tag_loc / my_notch_distance);
RF_theta = (tag_loc - closest_notch*my_notch_distance)/my_notch_distance*2*pi;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% we also need to compute the  gradient moment from the extra gradient area
% in the refocusing pulse (Gz_phaser).  This is what produces the phase accumulation
% the labeling location.  This is a function of the gradient and the
% location of the inversion pulse
% RF_theta = tag_loc * gambar * Gz_phaser * 1e-4 * vel*(0.25)^2
%
xtra_phase = gam * 1e-4 * A_xtra * tag_loc;  % in rads
RF_theta = -rem(xtra_phase, 2*pi) ;

RF_theta = RF_theta + phase_correction;
%RF_theta = 0;  % this line turns off the pre-phasing
%
%%
% calculate transmit frequency offset for desired tag location
% in order to achieve the inversion,
RF_dw = tag_loc * gambar * Gz_ss * 1e-4; % this is in kHz


%%
% integrate the pulses into the sequence - in Alsop's abstract: 500us
% pulses every 800 us

cnt = 0;
for t = 1/dt + 1: RF_spacing/dt : 400/dt

    if isTag
        doSign = 1;
    else
        doSign = (-1)^cnt;
    end

    % shift the frequency and add phase to the RF pulse
    FMx = cos(RF_dw  * 2*pi * [1:length(myHanning)]' * dt + cnt*RF_theta);
    FMy = sin(RF_dw  * 2*pi * [1:length(myHanning)]' * dt + cnt*RF_theta);
    cnt = cnt+1;

    myHanning_x = myHanning .* FMx;
    myHanning_y = myHanning .* FMy;

    % a train of RF pulses
    B1( t + t_ramp : ...
        t + t_ramp + length(myHanning) -1 ) = ...
        (myHanning_x  + i * myHanning_y) * doSign;

    % slice select gradient
    Gz( t : ...
        t + 4 * t_ramp + t_Gz_ss -1 ) = ...
        [trap_ss; trap_ss_r];
   
    
 end
% intorduce and artifactual gradient:
Gz = Gz + Gz_err;
%%

%B1(1100:1100+length(mySinc)-1) = mySinc;

% Throw in an ARTIFACTUAL gradient - B0 inhomogeneity
%Gz = Gz + 0.15;

% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz = (dx.*Gx + dy.*Gy + dz.*Gz) * 1e-4;

% these are the transverse fields - due to B1 pulse.
Bx = real(B1);
By = imag(B1);

%%

% Put the whole thing together and run the simulation
beff = [Bx By  Bz];
M = blochsim(Mi, beff, T1, T2, dt, nstep);
M = M(1:end-1,:);

%%
% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)
phi = gambar * sum(Bz * dt);

if ~doBatch

    %     figure
    %     subplot(411) , plot(Gx)
    %     subplot(412) , plot(Gy)
    %     subplot(413) , plot(Gz)
    %     subplot(413) , hold on, plot(1e5*real(B1) ,'r')
    %     subplot(414) , plot(real(B1)), hold on, plot(imag(B1),'r')

    figure

    subplot(311) , plot(M(:,1))
    subplot(312) , plot(M(:,2))
    subplot(313) , plot(M(:,3))


    figure
end
subplot(211),
plot(Gz(1:4/dt))
hold on
plot(1e5*real(B1(1:4/dt)) ,'r') , plot(1e5*imag(B1(1:4/dt)) ,'g')
hold off
subplot(212), plot(M(:,3))
axis([1 length(Gz) -1.1 1.1])
grid on



drawnow

if slomo
    figure
    range = -tag_loc/vel*dt
    for c=1.45e5:10:2.50e5
        subplot(311),compass(M(c,1), M(c,2))
        subplot(312), hold on; plot(c,Bx(c),'g*');
        subplot(313), hold on; plot(c,M(c,3),'g*');
        drawnow
    end
end
