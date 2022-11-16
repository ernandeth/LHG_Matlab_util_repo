function vel_sim_07
%
% Version uses the hyperbolic secant and tanget sweeps for robustness to
% off-resonance
%
% simulate an adiabatic pulse broken into chunks
% RF-bipolar-Rf-bipolar ...
% freqeuency of RF pulses sweeps from +deltaF to -deltaF
%
% This is the version with a sech pulse that's been broken up into chunks
% also, trying eddy current compensation.
%
% initial position, initial velocity, and acceleration
% units are cm and seconds (for now)

close all
% Display flags:
global showSequence showMag slomo showBeff doXlatebin doSinGrads off_resonance
showSequence = 1;
showMag = 1;
slomo = 0 ;
showBeff =0;
doXlatebin = 0;
doSinGrads = 1;
off_resonance =  0e-3 / 42576 ; % Hz gammabar is in kHz/T 

% default movement parmeters for spin
% start location, velocity and acceleration
movparms = [-1, 60, 0];

% Testing different pulse sequence parms:  make it shorter!!  This one
% works well

seqparms(1) = 0.2e-4;  % (default: 0.2e-4) B1 max (Tesla)  ;
seqparms(2) = 32 ; % N. of pulses in train ;
seqparms(3) = 2 ; % RF_duration ; % ms (doesn't matter)
seqparms(4) = 0.5 ;  % (default: 0.5) sweepWidth :  RF frequency sweep width (KHz)
seqparms(5) = 0 ;  % RF_phase correction
seqparms(6) = 2.4 ; % (2 ms default)Gz_duration; % duration in ms.(includes both Gz lobes)

seqparms(7) = 3  ; % Gz_amp ; % default = 3.5 Gauss / cm
seqparms(8) = 1;   % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(9) = 0;   % velocity target to invert


% sweep rate: kHz/ms =
sweepRate = seqparms(4) / (seqparms(2) * (seqparms(6)+seqparms(3)) );

% Things to try:
velocities =  linspace(-50,50,30);
velocities = 0;
alphas=[];

sweep_range = linspace(0.1, 3, 20);
sweep_range = seqparms(4);

starting_positions = [-10:10:10];

Gmax_arr = linspace(0,5,20);
Gmax_arr = seqparms(7);

B1max_arr = (1e-5 * linspace(5,0,5));
B1max_arr = seqparms(1);

Nseg_arr = round(linspace(1,48,20));
Nseg_arr = seqparms(7);

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
%for sweep = sweep_range
%      seqparms(4) = sweep;



% Try different number of segmentsf
%for Nseg =  Nseg_arr
%    seqparms(2) = Nseg;


% Try different B1max amplitudes
%for  B1max = B1max_arr
%    seqparms(1) = B1max

% Try diggerent gradient amplitudes
%for gm = Gmax_arr
%    seqparms(7) = gm;
    
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
%end

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
global showSequence showMag slomo showBeff doXlatebin doSinGrads off_resonance


% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;


% default target velocities for inversion
vel_target = 0 ; % in cm/s
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
% T2 = 100;
% T1 = 1100;

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
rfgap = 0.4/dt;

rfgap = 0.050/dt;

Gz_length = Gz_duration / dt + rfgap;
Gzseg = zeros( Gz_length, 1);

% now with ramps:
Gzseg = ones( Gz_length, 1);

tramp = floor(0.35/dt);
tflat = 0.5*(Gz_length - 4*tramp - rfgap);
 
Gzseg(1:4) = 0;
Gzseg(4+1   :   tramp+4) = linspace(0,1,tramp);
Gzseg(4+tramp + tflat +1 : 3*tramp+ tflat +4) = linspace(1,-1, 2*tramp);
Gzseg(4+3*tramp+ tflat +1 : 4+3*tramp+ 2* tflat  ) = -1 ;
Gzseg(4+3*tramp+ 2*tflat +1 : 4+4*tramp+ 2* tflat  ) =  linspace(-1,0, tramp) ;

Gzseg( end-rfgap+1 :end) = 0;


if doSinGrads
    Gzseg(:) = 0;
    lobelen = length(Gzseg) - rfgap ; %-2*rfgap;
    Gzseg(1 : lobelen) = sin(linspace(0,2*pi,lobelen));
end
%fprintf('ramp times are : %f ', tramp*dt);

% Gzseg(end/2+1:end) = -Gzseg(end/2:-1:1)*Gfac;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design a continuous 4 ms. sech pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 1;
k = 1;

RF_duration = 4.0 / Npulses;
sechlen = 4.0 / dt;
B1seglen = sechlen/Npulses;

% amplitude of the sech pulse
sech_amp = sech(beta * linspace(-pi , pi, sechlen));
%sech_amp = cos(linspace(-pi/2 , pi/2, B1seglen*Npulses));
%sech_amp = sech_amp / norm(sech_amp*dt);
sech_amp = sech_amp / max(sech_amp);

% NOTE:  here B1max is not flip angle at all, but B1 max in Tesla
sech_amp =  B1max * sech_amp;

% frequency sweep of the sech pulse
sech_freq = -tanh(k*beta * linspace(-pi , pi, sechlen)); % works well

% sech_freq = sin(linspace(-pi/2 , pi/2, B1seglen*Npulses));
% sech_freq = tan(0.4*linspace(-pi , pi, B1seglen*Npulses));

sech_freq = sech_freq * 2/( max(sech_freq)-min(sech_freq));
sech_freq = sweepWidth * sech_freq;

RF_dw = sech_freq;

% time vector for the sech pulses
sech_tvec = RF_duration * Npulses * linspace(-0.5, 0.5, sechlen);
mysech = sech_amp .* exp(i * 2 * pi * sech_freq .* sech_tvec);

% or  ...
sech_phase = cumsum(sech_freq)*dt*2*pi;
mysech = sech_amp .* exp(i .* sech_phase);

%fprintf('\nDuration of the Original sech pulse: %f ms', Npulses*RF_duration);
%fprintf('\rMaximum B1 amp in my_sech %f Tesla', max(abs(sech_amp)));


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
Gzpad = zeros(B1seglen,1);

Gzmodule = [Gzpad ; Gzseg];

mod_length = length(Gzmodule);
mod_duration = RF_duration + Gz_duration;



RF_phases = zeros(Npulses,1);
start_phi = 0;
RF_segment_type = 6;

for n = 1: Npulses
    % version 6:  a regular sech that has been broken into chunks
    % this means that the frequency varies throughout each the pulse.
    B1seglen=floor(B1seglen);
    FM = mysech( (n-1)*B1seglen + [1:B1seglen] )';
    
    
    % add phase for velocity selection
  %  Mxy_phase = - gam * cumsum( zpos_target .* Gzmodule_tmp) * dt   ;
  %  RF_phase = Mxy_phase(end) ;
  %  FM = FM * exp(i*n* RF_phase  ) ;
    
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
    if n<Npulses  % don't put a gradient after the last segment
        Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule ; % * (-1)^(round((n+3)/Npulses)) ;
    end
    
end

% calculate the position if it moves with the target velocity and
% acceleration during the pulse train 
ttmp  = dt*[0:length(Gz)-1]';
zpos_target = pos0 + vel_target*ttmp + 0.5*accel_target*(ttmp.^2);

% add corresponding phase for velocity selection
vel_phase = -gam * cumsum( zpos_target .* Gz) * dt   ;
B1 = B1 .* exp(i*vel_phase  ) ;
  



% save to a file for the scanner...
save vsi_train_freqs.txt  RF_dw -ascii
save sim_parms.txt seqparms -ascii
save vsi_train_phases.txt  RF_phases -ascii






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
% off_resonance =  50e-3 /gambar;  % gambar is in kHz/T
Bz = Bz + off_resonance ; % introduce off-resonance

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
    ylabel('Amplitude (T)')
    
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
    area(t, Gz * 1e4)
    axis ([0 Npulses*modlength*dt -5 5])
   % hold on
   % plot( t, 10*abs(B1),'r--' );

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

%drawnow


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
tmp(isnan(tmp)) = 0;
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
