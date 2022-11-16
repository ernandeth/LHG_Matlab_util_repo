function vel_sim_06
%
% Version uses the hyperbolic secant and tanget sweeps for robustness to
% off-resonance
%
% simulate an adiabatic pulse broken into chunks
% RF-bipolar-Rf-bipolar ...
% freqeuency of RF pulses sweeps from +deltaF to -deltaF
%
% This is the version with the sinc excitation and rect acceleration
% profile

% initial position, initial velocity, and acceleration
% units are cm and seconds (for now)

close all
% Display flags:
global showSequence showMag slomo showBeff doXlatebin
showSequence = 1;
showMag = 0;
slomo = 0 ;
showBeff =0;
doXlatebin = 1;

% default movement parmeters for spin
% start location, velocity and acceleration
movparms = [-1, 60, 0];

% Testing different pulse sequence parms:  make it shorter!!  This one
% works well

seqparms(1) = 0.2e-4;  % B1 max (Tesla)  ;
seqparms(2) = 20 ; % N. of pulses in train ;
seqparms(3) = 2 ; % RF_duration ; % ms (doesn't matter)
seqparms(4) = 0.5 ;  % sweepWidth :  RF frequency sweep width (KHz)
seqparms(5) = 0 ;  % RF_phase correction
seqparms(6) = 2.0 ; % Gz_duration; % duration in ms.(includes both Gz lobes)

seqparms(7) = 2.8  ; % Gz_amp ; % Gauss / cm
seqparms(7) = 3.5  ; % Gz_amp ; % Gauss / cm
seqparms(8) = 1;   % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(9) = 10;   % velocity target to invert


% sweep rate: kHz/ms =
sweepRate = seqparms(4) / (seqparms(2) * (seqparms(6)+seqparms(3)) );

% Things to try:
velocities =  linspace(-50,50,30);
alphas=[];

sweep_range = linspace(0.1, 3, 20);
sweep_range = 0.5;

starting_positions = [-10:10:10];

Gmax_arr = linspace(0,5,20);
Gmax_arr =3.5;

B1max_arr = (1e-5 * linspace(3,0,20));
B1max_arr = (0.2e-4);

Nseg_arr = round(linspace(1,48,20));
Nseg_arr = 32;

init_pos = [-15:5:15];
man_phase_corr = linspace(0, pi/10, 50);



% Try several starting positions
% for pos = starting_positions
%        movparms(1) = pos;

% TRy different starting positions of the spin
% for pos = init_pos
%     movparms(1) = pos;


% Try different Manual Phase corrections;
% for mpc = man_phase_corr;
%     seqparms(5) = mpc;

% Try different flip durations
% for  RF_dur=[0.5]
%    seqparms(3) = RF_dur;




%velocities = 0;

my_inversions = [];
my_controls = [];

% Try different sweep widths
% for sweep = sweep_range
%      seqparms(4) = sweep;



% Try different number of segments
%for Nseg =  Nseg_arr
%    seqparms(2) = Nseg;
    

% Try different B1max amplitudes
%for  B1max = B1max_arr
%    seqparms(1) = B1max

% Try diggerent gradient amplitudes
for gm = Gmax_arr
    seqparms(7) = gm;
    
    inversion = zeros(size(velocities));
    control = zeros(size(velocities));
    
    for n = 1:length(velocities)
        
        movparms(2) = velocities(n);
        %movparms(3) = -20;
        
        % EXecute the simulation here
        M = vel_sim(seqparms, movparms, 0, 0);
        doXlatebin = 0;
        
        % the z-magnetization at the end of the whole interval.
        inversion(n) = M(end,3);
        
        %Mc = vel_sim(seqparms, movparms, 1, 0);
        %control(n) = Mc(end,3);
    end
    
    my_inversions = [my_inversions ; inversion];
    my_controls = [my_controls ; control];
    
    if length(velocities) > 1
        
        figure(18)
        hold on
        plot(velocities, inversion)
        title('Final Mz as function of velocity.')
        axis([min(velocities)  max(velocities) -1 1]);
        xlabel('Velocity (cm/s)')
        ylabel('M_z')
        grid on
    end
    
    alphas = [alphas; min(inversion)];
    drawnow
end

%% Just making a nice figure to illustrate
% figure(84)
% subplot(313)
% plot(velocities, inversion)
% title('Velocity Profile of VSAI pulses.')
% xlabel('Velocity (cm/s)')
% ylabel('M_z')
% grid on


%movparms(2) = 30;

%M = vel_sim(seqparms, movparms, 0, 1);

if length(B1max_arr) > 1 &&  length(velocities) > 1
    figure(86)
    imagesc(my_inversions);
    xlabel('Velocities (cm/s)')
    
    set(gca,'XTick', [1 15.5 30]) ;
    set(gca,'YTick', [1 10.5 20]);
    set(gca,'XTickLabel', [-50 0 50]);
    set(gca,'YTickLabel', [300 150 0]);
    title('Longitudinal Magnetization');
    colormap jet
    colorbar
    caxis([ -1 1])

    xlabel('Velocities (cm/s)')

    ylabel('B_1max (mGauss)')
    title('M_z after VSAI pulse')
    colormap jet
    colorbar
 
end

if length(sweep_range) > 1 &&  length(velocities) > 1
    figure(87)
    imagesc(my_inversions);
    
    set(gca,'XTick', [1 15.5 30]) ;
    set(gca,'YTick', [1 10.5 20]);
    set(gca,'XTickLabel', [-50 0 50]);
    set(gca,'YTickLabel', [100 1550 3000]);
    title('Longitudinal Magnetization');
    colormap jet
    colorbar
    caxis([ -1 1])

    xlabel('Velocities (cm/s)')
    ylabel('Sweep Width (Hz)')
    
    title('M_z after VSAI pulse')
    colormap jet
    colorbar
 
end

if length(Nseg_arr) > 1 &&  length(velocities) > 1
    figure(88)
    imagesc(my_inversions);
    
    set(gca,'XTick', [1 15.5 30]) ;
    set(gca,'YTick', [1  20]);
    set(gca,'XTickLabel', [-50 0 50]);
    set(gca,'YTickLabel', [1  48]);
    title('Longitudinal Magnetization');
    colormap jet
    colorbar
    caxis([ -1 1])

    xlabel('Velocities (cm/s)')
    ylabel('N. of segments')
    
    title('M_z after VSAI pulse')
    colormap jet
    colorbar
 
end

if length(Gmax_arr) > 1 &&  length(velocities) > 1
    figure(88)
    imagesc(my_inversions);
    
    set(gca,'XTick', [1 15.5 30]) ;
    set(gca,'YTick', [1 10.5 20]);
    set(gca,'XTickLabel', [-50 0 50]);
    set(gca,'YTickLabel', [0 2.5 5]);
    title('Longitudinal Magnetization');
    colormap jet
    colorbar
    caxis([ -1 1])

    xlabel('Velocities (cm/s)')
    ylabel('G_{max} (G/cm)')
    
    title('M_z after VSAI pulse')
    colormap jet
    colorbar
 
end

    % figure(89)
% plot(sweep_range , alphas)
% plot(gmax_range , alphas)
% plot(man_phase_corr, alphas);


% movparms(2) = 50;
%
% M = vel_sim(seqparms, movparms, 0);

save velsim06_targetvel.mat

return



function M = vel_sim(seqparms, movparms, isControl , makeNiceFig);
%%
% display flags
global showSequence showMag slomo showBeff doXlatebin


% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;


% default target velocities for inversion
vel_target = 10 ; % in cm/s
accel_target = 0;  % in cm/s2

% default movement parameters
accel = 0 ; % in cm/s^2
vel = 55  ; % in cm/s
pos0 = -5 ;  % cm


% default pulse sequence parameters
Npulses = 12;
RF_duration = 0.5 ; % ms
sweepWidth = 0;     % RF frequency sweep width
RF_phase = 0;  % RF phase
Gz_duration = 2; % duration in ms.(includes both lobes)
Gz_amp = 1 ; % Gauss / cm
Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm
Gfac = 1.2;

%% overrride defaults if the user gives you parms
if nargin>0
    
    accel = movparms(3) ;
    vel = movparms(2) ;
    pos0 = movparms(1);
    
    
    % RF pulse parameters
    B1max = seqparms(1);
    Npulses = seqparms(2);
    RF_duration = seqparms(3) ; % ms
    
    sweepWidth = seqparms(4);     % RF frequency sweep width
    MPC = seqparms(5) ; % RF phase
    
    % gradient pulses parameters
    Gz_duration = seqparms(6); % duration in ms.(includes all fours Gz lobes)
    Gz_amp = seqparms(7) ; % Gauss / cm
    Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm
    
    Gfac = seqparms(8);  % fraction of the gradient
    
    vel_target = seqparms(9);   % Target velocity to invert
    
    seqparms = seqparms';
else
    
    % default pulse sequence parameters
    B1max = 35e-7;  % Tesla = 35 mG
    Npulses = 12;
    RF_duration = 0.5 ; % ms
    sweepWidth = 0;     % RF frequency sweep width
    RF_phase = 0;  % RF phase
    Gz_duration = 3*4; % duration in ms.(includes all fours Gz lobes)
    Gz_amp = 1 ; % Gauss / cm
    Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm
    Gfac = 1.2;
end
%%

% units: change from seconds to ms
accel = accel * 1e-6;
accel_target = accel_target * 1e-6;
vel = vel * 1e-3;
vel_target = vel_target * 1e-3;

% The usual constants
gambar = 42576;     % gamma in kHz/T
gam = gambar*2*pi;  % radians / (ms *Tesla)
mygam = gam / 1000; % radians / (s*Tesla)

% blood
T1 = 1660;  %ms
T2 = 175;   %ms

% brain
%T2 = 100;
%T1 = 1100;

B0 = 3;    % keep things in Tesla

% total duration of the simulation interval (in ms)
duration = 100;
NSTEPS = duration / dt;  % duration ms @  10 us steps.
t = linspace(0,duration, NSTEPS)';

Mi = [0, 0, 1]; % initial condition: tipped into x axis.

%%
% initialize gradients in Tesla/cm.
Gx = zeros(NSTEPS,1);
Gy = zeros(NSTEPS,1);
Gz = zeros(NSTEPS,1);

% initialize the B vector (tesla)
Bx = zeros(NSTEPS,1);
By = zeros(NSTEPS,1);
Bz = zeros(NSTEPS,1);

% initialize B1 channel
B1 = zeros(NSTEPS,1);

% Constructing the Bipolar Gradient pulses:
% Gradient is a 1  -1
Gz_length = Gz_duration / dt;
Gzseg = zeros( Gz_length, 1);

% adding a little extra time to prevent overlaps. (20 us)
%Gzseg(1:end/2) = -Gfac;
Gzseg(1+4 : end/2) = -Gfac; 
Gzseg = Gzseg * Gz_amp;


% now with ramps:
Gz_length = Gz_duration / dt;
Gzseg = ones( Gz_length, 1);
tramp = floor(Gz_length/4);

% adding a little extra time to prevent overlaps. (20 us)
%Gzseg(1:tramp) = linspace(0,1,tramp);
Gzseg(1+4 : tramp+4) = linspace(0,1,tramp);
Gzseg(1:4) = 0; 

%Gzseg(tramp+1 : 2*tramp) = linspace(1,0,tramp);
Gzseg(tramp +1 +4 : 2*tramp+4) = linspace(1,0,tramp);

if Gz_duration >1
    Gzseg(:) = 1;
    tramp = 0.4/dt;

    % adding a little extra time to prevent overlaps. (20 us)
    %Gzseg(1:tramp) = linspace(0,1,tramp);
    Gzseg(4+1 : tramp+4) = linspace(0,1,tramp);
    Gzseg(1:4) = 0;
    
    Gzseg(end/2 - tramp +1: end/2) = linspace(1,0, tramp);
end


Gzseg(end/2+1:end) = -Gzseg(1:end/2)*Gfac;

Gzseg = Gzseg * Gz_amp;
Gzpad = zeros(RF_duration/dt,1);
Gzmodule = [Gzpad ; Gzseg];
Gzmodule_tmp = Gzmodule;

% Now calculate the position of the particle over time:
% a time vector for the duration of one module
ttmp  = dt*[0:length(Gzmodule)-1]';

% calculate the position if it moves with the target velocity and acceleration:
zpos_target = pos0 + vel_target*ttmp + 0.5*accel_target*(ttmp.^2);


% calculate the phase gained by moving spins during one module
% integrate the position(t) * G(t)
%
% gam: rad/ms/Tesla
% Gz: T/cm
% dt : ms
Mxy_phase = - gam * cumsum( zpos_target .* Gzmodule) * dt   ;
RF_phase = Mxy_phase(end) ;

% add the manual phase correction
RF_phase = RF_phase + MPC ;

% control case:  add  +pi to the phase increment;
% this will cause every other RF pulse to go in the opposite direction
if isControl
    % RF_phase = Mxy_phase(end) + pi;
    % RF_phase = Mxy_phase(end) - MPC;
    % RF_phase  = RF_phase + pi;
    % RF_phase = -RF_phase;
    % Gzseg = -Gzseg;
    % Gzseg(:) = 0;
    % RF_phase = 0;
    Gzseg = -Gzseg;
    %RF_phase = -RF_phase;
end

% alternative
% if isControl
%     RF_phase = -RF_phase;
% end

%%
% constructing the RF pulse SEGMENTS with unit area:

% define a 0.5 ms Hanning pulse of unit area
segment_han = hanning (RF_duration/ dt);
segment_han = segment_han / sum(segment_han*dt);
hlen = length(segment_han);

% scale it to achieve max. flip angle of the train
% the area required for flip angle is = flip / gamma
% amplitude is in Tesla units
% B1max must be in rads
% gam is in rad/ms/Tesla

B1area = B1max  / gam;  % area of B1 in (Tesla * ms)

B1amp_hanning = B1area * segment_han;  % units are Tesla

segment_rect = ones(size(segment_han));
segment_rect = segment_rect / sum(segment_rect*dt);
B1amp_rect = B1area * segment_rect;  % units are Tesla

segment_sr = softrect(length(segment_han))';
segment_sr = segment_sr / sum(segment_sr*dt);
B1amp_sr = B1area * segment_sr;  % units are Tesla


% default amplitude modulation
B1amp = B1amp_hanning;

% weights and frequencies of sech pulse segments in the train

% amplitudes for each pulse in train
beta = 1;  % works well
beta = 1;
k = 0.5;

weights = sech(beta * linspace(-pi , pi, Npulses));
weights = weights - min(weights);
weights = weights / (max(weights) - min(weights));
%weights(:) =1;

% frequency for each of the pulses
RF_dw =  tanh(k *beta * linspace(-pi , pi, Npulses)); % works well
%RF_dw = sweepWidth * tanh(0.3*beta * linspace(-pi , pi, Npulses)); % works well
RF_dw = sweepWidth * RF_dw * 2/( max(RF_dw)-min(RF_dw));
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ... or a continuous 4 ms. sech pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 1;
k = 1;

RF_duration = 4.0 / Npulses;
sechlen = 4.0 / dt;
hlen = sechlen/Npulses;

% amplitude of the sech pulse
sech_amp = sech(beta * linspace(-pi , pi, sechlen));
%sech_amp = cos(linspace(-pi/2 , pi/2, hlen*Npulses));
%sech_amp = sech_amp / norm(sech_amp*dt);
sech_amp = sech_amp / max(sech_amp);

% NOTE:  here B1max is not flip angle at all, but B1 max in Tesla
sech_amp =  B1max * sech_amp;   

% frequency sweep of the sech pulse
sech_freq = -tanh(k*beta * linspace(-pi , pi, sechlen)); % works well

%sech_freq = sin(linspace(-pi/2 , pi/2, hlen*Npulses));
% sech_freq = tan(0.4*linspace(-pi , pi, hlen*Npulses));

sech_freq = sech_freq * 2/( max(sech_freq)-min(sech_freq));
sech_freq = sweepWidth * sech_freq;

RF_dw = sech_freq;

% time vector for the sech pulses
sech_tvec = RF_duration * Npulses * linspace(-0.5, 0.5, sechlen);
mysech = sech_amp .* exp(i * 2 * pi * sech_freq .* sech_tvec);

% or  ...
% sech_phase = cumsum(sech_freq)*dt*2*pi;
% mysech = sech_amp .* exp(i .* sech_phase);

%fprintf('\nDuration of the Original sech pulse: %f ms', Npulses*RF_duration);
%fprintf('\rMaximum B1 amp in my_sech %f Tesla', max(abs(sech_amp)));

%% or 8 ms BIR pulse
%{
Npoints = 4.0/dt;  % NUmber of points in original 90 pulse
Tseg = 4e-3;  % seconds  : duration of the original 90 pulse
tt = linspace(0,Tseg,Npoints);
zeta = 15;
kappa = atan(60); % 1.5541;
tankappa = 60;
wmax = 0.389e3 *2*pi;  % rad/sec  
wmax = 2*pi*1e3*sweepWidth/2;
amp1 = tanh(zeta * (1 - (tt./Tseg))); 
amp2 = tanh(zeta * (tt./Tseg)); 
sech_amp = [amp2 amp1];
sech_amp = sech_amp / sum(sech_amp*dt);
sech_amp =  B1area * sech_amp;

phase1 = -wmax*Tseg*log( abs(cos(kappa*tt/Tseg)) / (kappa*tankappa) );
phase2 = -wmax*Tseg*log( abs(cos(kappa*(tt/Tseg -1))) / (kappa*tankappa) );
sech_phase = [phase2 phase1];

mysech = sech_amp .* exp(i * sech_phase);
%}
%%

%{
    subplot(311)
    plot(abs(mysech)), title('Amplitude')
    axis tight
    hold on
    subplot(312)
    plot(sech_freq), title('Frequency')
        axis tight
    hold on
    subplot(313)
    plot(angle(mysech)), title('Phase')
    hold on
        axis tight
%}



% For stitching the elements together
% need to create some zero padding
B1pad = zeros(size(Gzseg));
Gzpad = zeros(hlen,1);

B1module = [B1amp ; B1pad];
Gzmodule = [Gzpad ; Gzseg];
mod_length = length(B1module);
mod_duration = RF_duration + Gz_duration;



RF_phases = zeros(Npulses,1);
start_phi = 0;
RF_segment_type = 6;

for n = 1: Npulses
    % step size is in usec. frequency is in kHz.
    % RF_duration is in ms.
    
    % generate each RF segment
    switch RF_segment_type
        
        case 1
            % version 1: constant FM at each pulse, hanning AM
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = weights(n)*B1amp_hanning;
            
        case 2
            % version 2:  hanning FM and hanning AM
            % pulse segments from 0 to RF_dw and back to 0
            RF_dw2 = RF_dw(n) * hanning(length(segment_han));  % kHz
            B1amp = weights(n)*B1amp_hanning;
            
        case 3
            % version 3:  constant FM and constant AM
            % pulse segments from 0 to RF_dw and back to 0
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = weights(n)*B1amp_rect;
            
        case 4
            % version 4: constant FM at each pulse, soft rect AM
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = weights(n)*B1amp_sr;
            
        case 5
            % version 5: soft rect  FM at each pulse, soft rect AM
            RF_dw2 = RF_dw(n) * segment_sr / max(segment_sr);  % kHz
            B1amp = weights(n)*B1amp_sr;
            
        case 6
            % version 6:  a regular sech that has been broken into chunks
            % this means that the frequency varies throughout each the pulse.
            hlen=floor(hlen);
            FM = mysech( (n-1)*hlen + [1:hlen] )';
            
    end
    
    
    
    
    % add phase for stitching together RF pulses without jumps
    if RF_segment_type ~= 6
        % B1 amp is the amplitude modulation in Tesla
        % RF_dw2 is the frequency modulation in kHz
        
        FM =  B1amp .* exp(i *2*pi*  RF_dw2 .* linspace(0, RF_duration, length(segment_han) )');
        
        % phase stitching here
        RF_phases(n) = start_phi;
        FM = FM * exp(i * (-angle(FM(1)) + start_phi)) ;
        start_phi = angle(FM(end)) ;
        
    end
    
    % add phase for velocity selection
    Mxy_phase = - gam * cumsum( zpos_target .* Gzmodule_tmp) * dt   ;
    RF_phase = Mxy_phase(end) ;
    FM = FM * exp(i*n* RF_phase  ) ;
    %Gzmodule_tmp = -Gzmodule_tmp;
    
    % break it down into Real and Imaginary
    FMx = real(FM);
    FMy = imag(FM);
    
    
    % RF pulses can be weighted independently
    % units are Tesla
    B1seg = (FMx  + i * FMy);
    
    % update the excitation module
    B1module = [B1seg ; B1pad];
    Gzmodule = [Gzpad ; Gzseg];
    modlength = length(B1module);
    
    B1( (n-1)*modlength + 1 : n*modlength) = B1module;
    % 12/21/15:  note the alternating gz modules
    Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule; % * (-1)^(n-1) ; 
    
    
    % save to a file for the scanner...
    save vsi_train_amps.txt weights -ascii
    save vsi_train_freqs.txt  RF_dw -ascii
    save sim_parms.txt seqparms -ascii
    save vsi_train_phases.txt  RF_phases -ascii
    
    
end



seqEnd = n*modlength;
pulseDuration = dt * Npulses * modlength;

% fprintf('\nDuration of the Final pulse: %f ms', Npulses*modlength*dt);

if doXlatebin
     genScannerPulses(B1(1:seqEnd), Gz(1:seqEnd), dt );
end

% alternative control
% if isControl
%     Gz = -Gz;
%     Gz(:) = 0;
% end


% calculate the actual position of particle in one module
zpos = pos0 + vel*t + accel*(t.^2);

% position of particle in cm over the whole period
dx = 0;
dy = 0;
dz = zpos;



% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz = (dx.*Gx + dy.*Gy + dz.*Gz) ;


% off resonance!!
% Bz = Bz .* (1 + 0.05 * randn(size(Bz)));
% Bz = Bz .* (1 + 0.05*(2/3) * randn(size(Bz)));
% Bz = Bz  + 1e-5 * randn(size(Bz));
% Bz = Bz + 1e-5 ; % introduce off resonance

% B1 error  !!
% B1error = 1 + 0*randn(1);

% these are the transverse fields  due to B1 pulse (in Tesla )
Bx = real(B1) ; %* B1error;
By = imag(B1) ; %* B1error;

%%

% Put the whole thing together and run the simulation
beff = [Bx By  Bz];
M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
%M = M(1:end-1,:);

%%
% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)
phi = gambar * sum(Bz * dt);

if showBeff
    
    
    % % The usual constants
    % gambar = 42576;     % gamma in kHz/T
    % gam = gambar*2*pi;  % radians / (ms *Tesla)
    % mygam = gam / 1000; % radians / (s*Tesla)
    
    B_eff = sech_amp + i*sech_freq/gambar;
    mn = -max(abs(B1));
    mx = max(abs(B1));
    
    prec_rate = gam*abs(B_eff);   % in rad/ms
    prec_rate = prec_rate(1:end-1);
    
    tip_rate = abs(diff(angle(B_eff))) / dt;
    adiabticity = prec_rate ./ tip_rate;
    
    figure(132)
    subplot(221)
    plot(sech_amp);
    axis tight; axis square;
    title('Amplitude')
    
    subplot(222)
    plot(sech_freq)
    axis tight; axis square;
    title('Frequency')
    
    subplot(223)
    plot(B_eff); axis square;
    %axis([mn mx mn mx])
    title('B_{eff} Trajectory')
    
    subplot(224)
    plot(adiabticity)
    axis tight; axis square;
    title('Adiabaticity Trajectory')
    
end


% Make nice plots here
if showMag
    
    figure(3)
    set(gcf,'Name', 'Magnetization (x,y,z)')
    set(gcf,'Position',[ 600 100 500 600 ])
    timevec = [0:(length(M)-1)]*dt;
    
    subplot(311) , plot(timevec, M(:,1))
    ylabel('M_x')
    axis([0 pulseDuration -1 1]);
    subplot(312) , plot(timevec, M(:,2))
    axis([0 pulseDuration -1 1]);
    ylabel('M_y')
    subplot(313) , plot(timevec, M(:,3))
    ylabel('M_z')
    axis([0 pulseDuration -1 1]);
    xlabel ('Time (ms)')
    %subplot(414), plot(timevec(1:length(zpos)),zpos);
    
    grid on
    %
    %     figure(4)
    %     plot(t,dz)
    %     title('Position of particle in time')
    %
%     if abs(vel_target - vel) < 0.002
%         keyboard
%     end
        
end


if showSequence
    
    figure(1)
    set(gcf,'Name', 'Pulse sequence ')
    
    subplot(311)
    plot(t, abs(B1)) ;    
    title('B_1 Amplitude')
    axis ([0 Npulses*modlength*dt 0 (max(abs(B1))*1.1+eps) ])
    ylabel('Amplitude (mG)')
    
    %subplot(412)
    %plot(t, RF_dw) ;    title('B_1 Frequency')
    %axis tight
    %ylabel('Frequency (KHz)')
    
    subplot(312)
    plot( t, angle(B1) );
    title('Phase')
    axis ([0 Npulses*modlength*dt -pi pi])
    ylabel('Phase (rad.)')
    
    subplot(313)
    area(t, Gz)
    axis ([0 Npulses*modlength*dt -5e-4 5e-4])
    title('Grad')
    ylabel('G_z (G/cm)')
    xlabel('Time (ms)')
    
    
end

if makeNiceFig
    figure(84)
    
    subplot(311)
    area(t, Gz)
    hold on
    plot(t, 100*real(B1) ,'r') ,
    plot(t, 100*imag(B1) ,'g')
    hold off
    axis([0 5 -3e-4 3e-4])
    legend('Gz',  'B1x', 'B1y')
    xlabel('Time (ms)')
    title('VSAI  pulse train ')
    
    
    timevec = [0:(length(M)-1)]*dt;
    subplot(312) , plot(timevec, M(:,3))
    axis([0 timevec(end) -1 1]);
    title('Evolution of M_z under VSAI pulses')
    xlabel('Time (ms)')
    ylabel('M_z')
    
end


if slomo
    figure(2)
    
    for c=1:100:length(M)
        clf
        line([0 M(c,1)], [0 M(c,2)], [0 M(c,3)])
        axis([-1 1 -1 1 -1 1])
        drawnow
    end
end

drawnow


M = M(end,:);


return

function result = softrect(Npts)
%%
r = ones(1,Npts);
h = hanning(Npts/2);
sr = conv(r,h);

Npts2 = length(sr);
inds = round(linspace(1,Npts2,Npts));
result = sr(inds);
%%
return


function test_sechs
%%

gambar = 42576;     % gamma in kHz/T
sweepWidth = 1;
beta = 1;  % works well
beta = 2;

subplot(311), hold off
subplot(312), hold off
subplot(313), hold off

Npulses = 16;
beta = 1;

for k = 0.1:0.1:2
    weights = sech(beta * linspace(-pi , pi, Npulses));
    weights = weights - min(weights);
    weights = 5e-5 * weights / (max(weights) - min(weights));
    
    
    % frequency sweep of the pulse train envelope
    RF_dw = sweepWidth * tanh(k*beta * linspace(-pi , pi, Npulses)); % works well
    RF_dw = sweepWidth * RF_dw * 2/( max(RF_dw)-min(RF_dw));
    
    B_eff = complex(gambar*weights, RF_dw);
    
    prec_rate = 2*pi*gambar*abs(B_eff);   % in kHz
    prec_rate = prec_rate(1:end-1);
    tip_rate = 2*pi*diff(angle(B_eff)) / (0.5);
    adiabty = prec_rate ./ tip_rate;
    
    min(adiabty)
    
    subplot(311)
    plot(weights)
    hold on
    
    subplot(312)
    plot(RF_dw)
    hold on
    
    subplot(313)
    plot(adiabty)
    hold on
    
    pause(0.1)
    
    
end

return


function genScannerPulses(rf1, G, dt)
%
% this function makes the files for the GE scanner to use
% it calls xlatebin
%
% function genScannerPulses(rf1, G, dt)
% dt is in msec.
% the output files are in 4 usec rampling period

SAMPRATE = 0.004;
NPOINTS = floor(length(rf1) * dt / SAMPRATE);
dnsample = floor(SAMPRATE/dt);

rf1 = rf1(1:dnsample:end);
G = G(1:dnsample:end);

dacmax = hex2dec('7ffe');

% Now write out scanner files:  magnitude
tmp = abs(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
figure(1)
set(gcf,'Name', 'Pulse sequence ')

subplot(312)
plot( tmp ) ;    title('RF mag')


filename = sprintf('myVSI_%d.rho.txt',NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

% first write it as intt16
filename = sprintf('myVSI_%d.rho.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.rho myVSI_%d.rho.bin', NPOINTS, NPOINTS);
eval(str)

% Now write out scanner files:   phase
tmp = angle(rf1);
tmp = tmp * dacmax / max(abs(tmp));
tmp = 2*round(tmp/2);

subplot(313)
plot( tmp );    title('RF phase')


filename = sprintf('myVSI_%d.theta.txt',NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);

% first write it as intt16
filename = sprintf('myVSI_%d.theta.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.theta myVSI_%d.theta.bin', NPOINTS, NPOINTS);
eval(str)




% Now write out scanner files:
tmp = G ;
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
subplot(311)
area( tmp)
title('Grad')

% write out a text file with the gradient waveform
filename = sprintf('myVSI_%d.grad.txt',NPOINTS);
fp=fopen(filename,'w');
for n=1:NPOINTS, fprintf(fp,'%d\n',tmp(n)); end
fclose(fp);


% first write it as intt16
filename = sprintf('myVSI_%d.grad.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.grad myVSI_%d.grad.bin', NPOINTS, NPOINTS);
eval(str)

fprintf('\n...Generated GE wavefiles for myVSI...\n')



return
