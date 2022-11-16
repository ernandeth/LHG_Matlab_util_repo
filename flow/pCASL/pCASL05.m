% simulation for pseudoCASL experiment based on blochsim

% this section is skipped to allow for iterations when the doBatch flag is
% turned on.  YOu do this when you are calling the script from a batch
% file, e.g. - see pCASL_batch01.m
%
% this version: try to make the inversion velocity selective by changing 
% the spacing between pulses.  ... note that the phase gain from velocity
% is very sensitive to t_Gz_rephaser
%

%%
doBatch = 1;

if ~doBatch
    tag_loc = 3;   
    flip_ang = 22.5;
    isTag = 1;
    vel = 20;
    slomo =0;
    eta = 0.25;   % works for all locs if eta = 0.25!!
    t_ramp = 0.02;  % ramp time in ms.
    Gz_err = 0;  % artefactual gradient.  accounting for susceptibility
    RF_spacing = 1.5 ;  % default = 0.8.  check below.
    vel_phase = 0;  % added phase to the RF pulses to change the velocity selectivity.
    
    % eta = 0.5 means balanced slice-select and refocuser : 
    % the refocuser area is half the ss area
    slomo = 0;
    close all
	Gz_ss = 0.1;  % slice select gradient in G/cm.
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
% RF_spacing = 0.80 ; % time between RF pulses in ms. default = 0.8
t_ramp = 0.02;  % ramp time in ms.


t_Gz_ss = 0.5 ; % ms
t_Gz_rephaser = 0.30 ; % ms

% areas of gradients:  G/cm * ms
A_Gz_ss = Gz_ss * t_Gz_ss;
A_Gz_rephaser = A_Gz_ss/2;  % this is what a typical balanced rephaser should be.  It gets ignored.

% GZ_ss2 is the refocuser gradient plus the xtra gradient for phase accrual
% Gz_xtra is the gradient imbalance used to gain some phase
% uses definition of eta from Wu et al 2007
A_Gz_ss2 = A_Gz_ss * (1-eta);
%A_xtra = A_Gz_ss2 - A_Gz_rephaser;
A_xtra = A_Gz_ss2 - A_Gz_ss;


% %%%%%%%%%%%%%%  screwing around
Gz_ss = 0;
%%%%%%%%%%%%%%%%

% just a couple of sanity checks
myfractionalmoment = A_Gz_ss2/ A_Gz_ss;
myeta = (A_Gz_ss2 - A_Gz_ss)/A_Gz_ss;

% calculate gradient strengths from areas:  G/cm
Gz_rephaser = A_Gz_rephaser/t_Gz_rephaser;
Gz_ss2 = A_Gz_ss2 / t_Gz_rephaser;
Gz_xtra = A_xtra / t_Gz_rephaser;

% unit conversions
t = linspace(0,500,nstep)'; % time vector in ms
t_Gz_ss = t_Gz_ss / dt; % ms to samples
t_Gz_rephaser = t_Gz_rephaser/ dt; % ms to  samples

flip_ang = flip_ang * pi /180;  % flip angle in radians
vel = vel * 1e-6 ; % cm/s scaled to cm/us
t_ramp = t_ramp /dt;

%%
Mi = [0, 0, -1]; % initial condition: tipped into x axis.

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

% define a 0.5 ms Hanning pulse and Normalize it so the area is 1.
myHanning = hanning (t_Gz_ss );
myHanning = myHanning / sum(myHanning);
% scale it to achieve flip angle - result is in Tesla
% the area required for flip_angle is  = flip / gamma
% and put it in Tesla units
myHanning = myHanning * flip_ang / (gam * dt);

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

% adding extra phase for velocity selectivity.
RF_theta = RF_theta + vel_phase;

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
    Gz( t + t_ramp : ...
        t + t_ramp + length(myHanning) -1 ) = ...
        Gz_ss;
   
    % rephaser for slice select
    Gz( t + t_ramp  + length(myHanning) + 2*t_ramp : ...
        t + t_ramp  + length(myHanning) + 2*t_ramp + t_Gz_rephaser -1 ) = ...
        -Gz_ss2;
	
	% one more "rephaser"
	Gz( t + 2*t_ramp  + length(myHanning) + 2*t_ramp + t_Gz_rephaser -1  : ...
        t + 2*t_ramp  + length(myHanning) + 2*t_ramp + 2*t_Gz_rephaser -1 ) = ...
        Gz_ss2;
	
    
 end
% intorduce and artifactual gradient:
Gz = Gz + Gz_err;
%%

%B1(1100:1100+length(mySinc)-1) = mySinc;

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
