% simulation for pseudoCASL experiment based on blochsim

% this section is skipped to allow for iterations when the doBatch flag is
% turned on.  YOu do this when you are calling the script from a batch
% file, e.g. - see pCASL_batch01.m
%
% this version is intended to simulate adjustments in eta to compensate for
% susceptibility gradients

pcolor='black';
%%
doBatch = 1;
showplot = 0;

if ~doBatch

    clear all
    %close all
    showplot = 1;
    pcolor='red';
    %hold on;
    %%%%%%%%%%%%% TAGGING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tag_loc = 0;   
    vel = 60;
    %%%%%%%%%%%%% TAGGING PULSE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flip_ang = 30; % in Degrees
    eta = 0.3; % (A+ - A-)/ A+
    isTag = 1;
    
    flip_ang = flip_ang * pi /180;
    RF_spacing = 0.80 ; % time between RF pulses in ms.
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    t_Gz_ss = 0.5 ; % ms
    t_Gz_rephaser = 0.200 ; % ms

    xtra_phase2 =0; %(rad)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%% GRADIENT ERROR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
    Gradient_error = 1; % 0: No gradient error | 1 : gradient error
    Gradient_error_shape = 1; %1: global  2: local
    Gz_err_amp = 0.1;     % amplitude of gradient error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%OFF-RESONANCE ERROR PARAMETERS%%%%%%%%%%%%%%%%
    off_res_error = 1;
    off_res_peak = 350; % in Hz
    off_res_width = 3; %cm - width of off-resonance pattern
    off_shape=1; % 1:rect , 2:hanning - shape of off-resonance pattern
    error_loc = 3 ; % 1:before; 2:first 1/4; 3:middle; 4:first 3/4; 5: last
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixerror = 0; % 0: No fixing| 1: eta update| 2: prephasing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    slomo = 0;
    %close all
end
%%
%%%%%%%%%%%%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
gambar = 42570;    % gamma/2pi in kHz/T
gam = gambar*2*pi;
T1 = 1664;%  ms @ 3T
T2 = 250;   %ms
B0 = 3;    % keep things in Tesla
dt = 1e-3; % 1 us expressed in millliseconds.  dt converts from samples to ms
nstep = 500 / dt;  % 500 ms @  1 us steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% areas of gradients:  G/cm * ms
A_Gz_ss = Gz_ss * t_Gz_ss;
A_Gz_ss2 = A_Gz_ss * (1-eta);
A_xtra = A_Gz_ss2 - A_Gz_ss;

% calculate gradient strengths from areas:  G/cm
Gz_ss2 = A_Gz_ss2 / t_Gz_rephaser;


%% Hesam Stuff
del_eta=Gz_err_amp*RF_spacing/A_Gz_ss;
%=========== Phase Compensation Variables====
taw1 = t_Gz_ss; %ms
taw2 = t_Gz_rephaser; %ms
delta = RF_spacing; %ms
Ga = Gz_err_amp; %G/cm
spin_vel = vel; % cm/s
G1 = Gz_ss; % G/cm
G2 = Gz_ss2; % G/cm
%============================================

error_Phase1 = gam * 1e-4 * Ga * tag_loc * delta;
error_Phase2 = gam * 1e-4 * Ga * spin_vel *(delta*taw1 - 3*taw1^2/4 - taw2^2/2 - taw1*taw2/2 + taw2*delta + delta^2/2)*1e-3;
error_Phase = error_Phase1 + error_Phase2; % In rad
Gz_ss2;
del_G2 = error_Phase/(gam*1e-4*(taw2 * tag_loc + spin_vel*(3*taw1*taw2/2+taw2^2)*1e-3));

switch fixerror
    case 1
        Gz_ss2 = Gz_ss2 + del_G2;
        A_Gz_ss2fixed = Gz_ss2*t_Gz_rephaser;
        neweta = (A_Gz_ss-A_Gz_ss2fixed)/A_Gz_ss;

    case 2
        xtra_phase2 = 2*pi*off_res_peak*RF_spacing*1e-3/2;%(rad)
end

%%
% unit conversions
t = linspace(0,500,nstep)'; % time vector in ms
t_Gz_ss = t_Gz_ss / dt; % ms to samples
t_Gz_rephaser = t_Gz_rephaser/ dt; % ms to  samples

 % flip angle in radians
%vel = vel * 1e-6 ; % cm/s scaled to cm/us
t_ramp = t_ramp /dt;


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
Bz_error = zeros(nstep,1);
Gz_error =  zeros(nstep,1);


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
dz = dz' * vel * 1e-6; % cm/s scaled to cm/us

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
myHanning = hanning (0.5 / dt);
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
xtra_phase = gam * 1e-4 * A_xtra * -1*tag_loc + xtra_phase2;  % in rads
RF_theta = -rem(xtra_phase, 2*pi) ;

%RF_theta = 0;  % this line turns off the pre-phasing
%
%%
% calculate transmit frequency offset for desired tag location
% in order to achieve the inversion,
RF_dw = -tag_loc * gambar * Gz_ss * 1e-4; % this is in kHz


%%
% integrate the pulses into the sequence - in Alsop's abstract: 500us
% pulses every 800 us

cnt = 0;
for t = 1/dt + 1: RF_spacing/dt : 499/dt  % 499=500-1 .. 500 ms: simulation time

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
    
 end
% intorduce and artifactual gradient:
if Gradient_error == 1;
    switch Gradient_error_shape 
        case 1
           Gz = Gz + Gz_err_amp;
        case 2
           Gz_error(find(dz>(tag_loc - off_res_width/2) & dz<(tag_loc + off_res_width/2)))=Gz_err_amp;
           Gz = Gz + Gz_error;
    end
end
    
%%

%B1(1100:1100+length(mySinc)-1) = mySinc;

% Throw in an ARTIFACTUAL gradient - B0 inhomogeneity
%Gz = Gz + 0.15;

% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz_ideal = (dx.*Gx + dy.*Gy + dz.*Gz) * 1e-4; %in Tesla

Bz_error_peak = off_res_peak*1e-3/gambar; % in Tesla


if off_res_error ==1

    if off_shape == 1
        Bz_error(find(dz>(tag_loc - off_res_width/2) & dz<(tag_loc + off_res_width/2)))=Bz_error_peak;
    end

    if off_shape == 2

        switch error_loc
            case 1
                ll = find(dz>(tag_loc - off_res_width) & dz<tag_loc);
                Bz_error(find(dz>(tag_loc - off_res_width) & dz<tag_loc))=Bz_error_peak.*hanning(length(ll));
            case 2
                ll = find(dz>(tag_loc - 3*off_res_width/4) & dz<(tag_loc + off_res_width/4));
                Bz_error(find(dz>(tag_loc - 3*off_res_width/4) & dz<(tag_loc + off_res_width/4)))=Bz_error_peak.*hanning(length(ll));
            case 3
                ll=find(dz>(tag_loc - 2*off_res_width/4) & dz<(tag_loc + 2*off_res_width/4));
                Bz_error(find(dz>(tag_loc - 2*off_res_width/4) & dz<(tag_loc + 2*off_res_width/4)))=Bz_error_peak.*hanning(length(ll));
            case 4
                ll=find(dz>(tag_loc - 1*off_res_width/4) & dz<(tag_loc + 3*off_res_width/4));
                Bz_error(find(dz>(tag_loc - 1*off_res_width/4) & dz<(tag_loc + 3*off_res_width/4)))=Bz_error_peak.*hanning(length(ll));
            case 5
                ll=find(dz>(tag_loc) & dz<(tag_loc + off_res_width));
                Bz_error(find(dz>(tag_loc) & dz<(tag_loc + off_res_width)))=Bz_error_peak.*hanning(length(ll));

        end


    end

    Bz = Bz_ideal +Bz_error;
else
    Bz = Bz_ideal;

end
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

% % % if ~doBatch
% % % 
% % %     %     figure
% % %     %     subplot(411) , plot(Gx)
% % %     %     subplot(412) , plot(Gy)
% % %     %     subplot(413) , plot(Gz)
% % %     %     subplot(413) , hold on, plot(1e5*real(B1) ,'r')
% % %     %     subplot(414) , plot(real(B1)), hold on, plot(imag(B1),'r')
% % % 
% % %     figure
% % % 
% % %     subplot(311) , plot(M(:,1))
% % %     subplot(312) , plot(M(:,2))
% % %     subplot(313) , plot(M(:,3))
% % % 
% % % 
% % %     figure
% % % end
% % % subplot(211),
% % % plot(Gz(1:4/dt))
% % % hold on
% % % plot(1e5*real(B1(1:4/dt)) ,'r') , plot(1e5*imag(B1(1:4/dt)) ,'g')
% % % hold off
% % % subplot(212), plot(M(:,3))
% % % axis([1 length(Gz) -1.1 1.1])
% % % grid on
% % % 
% % % 
% % % figure;

if showplot
    plot(-1*dz(end-1:-1:1),M(:,3),pcolor)
    ylim([-1.1 1.1])
    xlabel(sprintf('Z location (cm)   ===   Min(Mz) = %2.2g    Mz(end) = %2.2g',min(M(:,3)), M(end,3)));
    grid on
    drawnow
end
