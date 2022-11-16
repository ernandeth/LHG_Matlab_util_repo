function vel_sim_05
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

% default movement parmeters for spin
% start location, velocity and acceleration
movparms = [-1, 60, 0];

% Testing different pulse sequence parms:  make it shorter!!  This one
% works well
seqparms(1) = deg2rad(80);  %flip_ang ;
seqparms(2) = 16 ; % N. of pulses in train ;
seqparms(2) = 25 ; % N. of pulses in train ;
seqparms(3) = 0.25; % RF_duration ; % ms
seqparms(4) = 0.5 ;  % sweepWidth :  RF frequency sweep width (KHz)
seqparms(5) = 0 ;  % RF_phase correction
seqparms(6) = 2.0 ; % Gz_duration; % duration in ms.(includes both Gz lobes)
seqparms(6) = 2.2 ; % Gz_duration; % duration in ms.(includes both Gz lobes)
seqparms(7) = 2.5  ; % Gz_amp ; % Gauss / cm
seqparms(8) = 1;   % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(9) = 0;   % velocity target to invert


% sweep rate: kHz/ms =
sweepRate = seqparms(4) / (seqparms(2) * (seqparms(6)+seqparms(3)) );

% Things to try:
velocities =  linspace(-50,50,40);
alphas=[];
sweep_range = [0.1:0.1:1];
starting_positions = [-10:10:10];
gmax_range = [0.5:0.5:3];

flips = deg2rad(linspace(120,1,10));

init_pos = [-15:5:15];
man_phase_corr = linspace(0, pi/10, 50);

% Try diggerent gradient amplitudes
%for gm = gmax_range
%    seqparms(7) = gm;


% Try several starting positions
% for pos = starting_positions
%        movparms(1) = pos;


%velocities = 30;

% Try diggerent flip durations
%for  RF_dur=[0.5]
%    seqparms(3) = RF_dur;

my_inversions = [];
my_controls = [];

% Try different sweep widths
%for sweep = sweep_range
%     seqparms(4) = sweep;



% TRy different starting positions of the spin
% for pos = init_pos
%     movparms(1) = pos;


% Try different Manual Phase corrections;
% for mpc = man_phase_corr;
%     seqparms(5) = mpc;

% Try diggerent flip angles
for  fa = flips
    seqparms(1) = fa;
    
    inversion = zeros(size(velocities));
    control = zeros(size(velocities));
    
    for n = 1:length(velocities)
        
        movparms(2) = velocities(n);
        %movparms(3) = -20;
        
        M = vel_sim(seqparms, movparms, 0, 0);
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
        % plot(velocities, control, 'g')
        % plot(velocities, control-inversion, 'r')
        % legend('label','control','subtraction')
        title('Final Mz as function of velocity.')
        axis([min(velocities)  max(velocities) -1 1]);
        xlabel('Velocity')
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


figure(85)
imagesc(my_inversions);
xlabel('Velocities (cm/s)')
ylabel('Start position')
ylabel('Flip Angle (degrees)')
set(gca,'XTick', 1:4:length(velocities)) ;
set(gca,'YTick', 1:3:length(flips));
set(gca,'XTickLabel', round(velocities(1:4:end)));
set(gca,'YTickLabel', round(rad2deg(flips(1:3:end))));
title('Longitudinal Magnetization');
colormap jet
colorbar

% figure(86)
% imagesc(my_inversions - my_controls);
% xlabel('Velocities')
% ylabel('Start position')
% ylabel('Flip Angles')
% title('Subtractions');
% colormap jet

% figure(89)
% plot(sweep_range , alphas)
% plot(gmax_range , alphas)
plot(rad2deg(flips) , alphas);  title('Adiabatic Inversion'); xlabel('Max Flip'); ylabel('Final Mz')
% plot(man_phase_corr, alphas);


% movparms(2) = 50;
%
% M = vel_sim(seqparms, movparms, 0);

return



function M = vel_sim(seqparms, movparms, isControl , makeNiceFig);
%%
showSequence = 0;
showMag = 0;
slomo = 0;
showBeff = 0;


% target velocities for inversion
vel_target = 10 ; % in cm/s
accel_target = 0;  % in cm/s2

% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;

% default movement parameters
accel = 0 ; % in cm/s^2
vel = 55  ; % in cm/s
pos0 = -5 ;  % cm


% default pulse sequence parameters
flip_ang = deg2rad(18);
flip_ang = deg2rad(15);
Npulses = 12;
RF_duration = 0.5 ; % ms
sweepWidth = 0;     % RF frequency sweep width
RF_phase = 0;  % RF phase
Gz_duration = 2; % duration in ms.(includes both lobes)
Gz_amp = 1 ; % Gauss / cm
Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm
Gfac = 1.2;

% overrride defaults if the user gives you parms
if nargin>0
    
    accel = movparms(3) ;
    vel = movparms(2) ;
    pos0 = movparms(1);
    
    
    % RF pulse parameters
    flip_ang = seqparms(1);
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
    flip_ang = deg2rad(18);
    flip_ang = deg2rad(15);
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

T1 = 1660;  %ms
T2 = 200;   %ms

B0 = 3;    % keep things in Tesla

% total duration of the simulation interval (in ms)
duration = 500;
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
Gzseg = ones( Gz_length, 1);
Gzseg(1:end/2) = -Gfac;
Gzseg = Gzseg * Gz_amp;


% now with ramps:
Gz_length = Gz_duration / dt;
Gzseg = ones( Gz_length, 1);
tramp = floor(Gz_length/4);
Gzseg(1:tramp) = linspace(0,1,tramp);
Gzseg(tramp+1:2*tramp) = linspace(1,0,tramp);

if Gz_duration >1
    Gzseg(:) = 1;
    tramp = 0.2/dt;
    Gzseg(1:tramp) = linspace(0,1,tramp);
    Gzseg(end/2 - tramp +1: end/2) = linspace(1,0, tramp);
end


Gzseg(end/2+1:end) = -Gzseg(1:end/2)*Gfac;

Gzseg = Gzseg * Gz_amp;


% constructing the RF pulse SEGMENTS with unit area:

% define a 0.5 ms Hanning pulse 
segment_han = hanning (RF_duration/ dt);
hlen = length(segment_han);

segment_han = segment_han / sum(segment_han*dt);
 
% scale it to achieve max. flip angle of the train
% the area required for flip_angle is = flip / gamma
% amplitude is in Tesla units
% flip_ang must be in rads
% gam is in rad/ms/Tesla

B1area = flip_ang  / gam;  % area of B1 in (Tesla * ms) 

B1amp_hanning = B1area * segment_han;  % units are Tesla

segment_rect = ones(size(segment_han));
segment_rect = segment_rect / sum(segment_rect*dt);
B1amp_rect = B1area * segment_rect;  % units are Tesla

segment_sr = softrect(length(segment_han))';
segment_sr = segment_sr / sum(segment_sr*dt);
B1amp_sr = B1area * segment_sr;  % units are Tesla


% default amplitude modulation
B1amp = B1amp_hanning;

% Putting RF and Gz together into a TRAIN of pulses:

% amplitudes for the pulse train envelope
beta = 1;  % works well
beta = 1;
weights = sech(beta * linspace(-pi , pi, Npulses));
weights = weights - min(weights);
weights = weights / (max(weights) - min(weights));
%weights(:) =1;

% frequency sweep of the pulse train envelope 
RF_dw = sweepWidth * tanh(0.3*beta * linspace(-pi , pi, Npulses)); % works well
RF_dw = sweepWidth * tanh(0.3*beta * linspace(-pi , pi, Npulses));
RF_dw = sweepWidth * RF_dw * 2/( max(RF_dw)-min(RF_dw));

% For stitching the elements together
% need to create some zero padding
B1pad = zeros(size(Gzseg));
Gzpad = zeros(size(segment_han));

B1module = [B1amp ; B1pad];
Gzmodule = [Gzpad ; Gzseg];
mod_length = length(B1module);
mod_duration = RF_duration + Gz_duration;


% Now calculate the position of the particle over time:
% a time vector for the duration of one module
ttmp  = linspace(0 , mod_duration, mod_length)';

% calculate the position if it moves with the target velocity and acceleration:
zpos_target = pos0 + vel_target*ttmp + 0.5*accel_target*(ttmp.^2);


% calculate the phase gained by moving spins during one module
% integrate the position(t) * G(t)
%
% gam: rad/ms/Tesla
% Gz: T/cm
% dt : ms
Mxy_phase = -gam * cumsum( zpos_target .* Gzmodule)*dt   ;
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

RF_phases = zeros(Npulses,1);
start_phi = 0;
RF_segment_type = 4;

for n = 1: Npulses
    % step size is in usec. frequency is in kHz.
    % RF_duration is in ms.
    
    % generate each RF segment
    switch RF_segment_type
        
        case 1
            % version 1: constant FM at each pulse, hanning AM
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = B1amp_hanning;
            
        case 2
            % version 2:  hanning FM and hanning AM
            % pulse segments from 0 to RF_dw and back to 0
            RF_dw2 = RF_dw(n) * hanning(length(segment_han));  % kHz
            B1amp = B1amp_hanning;
            
        case 3
            % version 3:  constant FM and constant AM
            % pulse segments from 0 to RF_dw and back to 0
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = B1amp_rect;
        
        case 4
            % version 4: constant FM at each pulse, soft rect AM
            RF_dw2 = RF_dw(n) * ones(size(segment_han));  % kHz
            B1amp = B1amp_sr;
        
        case 5
            % version 5: soft rect  FM at each pulse, soft rect AM
            RF_dw2 = RF_dw(n) * segment_sr / max(segment_sr);  % kHz
            B1amp = B1amp_sr;
    end
    
    
    
    % B1 amp is the amplitude modulation in Tesla
    % RF_dw2 is the frequency modulation in kHz
    FM =  B1amp .* exp(i *2*pi*  RF_dw2 .* linspace(0,RF_duration,length(segment_han) )');

    % add phase for stitching together RF pulses without jumps
    RF_phases(n) = start_phi;
    FM = FM * exp(i * start_phi) ; 
    start_phi = angle(FM(end)) ;

    
    % add phase for velocity selection
    FM = FM * exp(i*n* RF_phase) ;

    % break it down into Real and Imaginary
    FMx = real(FM);
    FMy = imag(FM);
    
       
    % RF pulses can be weighted independently
    % units are Tesla
    B1seg = weights(n)*(FMx  + i * FMy);
    
    % update the excitation module
    B1module = [B1seg ; B1pad];
    Gzmodule = [Gzpad ; Gzseg];
    modlength = length(B1module);
    
    B1( (n-1)*modlength + 1 : n*modlength) = B1module;
    Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule;
    
    
    % save to a file for the scanner...
    save vsi_train_amps.txt weights -ascii
    save vsi_train_freqs.txt  RF_dw -ascii
    save sim_parms.txt seqparms -ascii
    save vsi_train_phases.txt  RF_phases -ascii
        
    if showBeff
        % B_eff for each segment of the pulse
        % (remember gambar is in kHz/T, B1seg is in Tesla)
        B_eff = abs(B1seg) + i*RF_dw2 / gambar;
                
        prec_rate = 2*pi*gambar*abs(B_eff);   % in kHz
        prec_rate = prec_rate(1:end-1);
        tip_rate = 2*pi*diff(angle(B_eff)) / dt;
        adiabty = prec_rate ./ tip_rate;
        
        
        figure(311);
        
        subplot(231)
        plot(B1amp * weights(n)); title('|b1|')
        subplot(232)
        plot(RF_dw2); title('pulse freq'); axis([0 length(B1seg) -sweepWidth sweepWidth])
        subplot(233)
        plot(B_eff, '*') ; axis([-1 1 -1 1]*5e-5); axis square
        title('B effective (seg)'); xlabel('B1'), ylabel('delta W')
        
        subplot(234)
        plot(weights);
        title('B1 weights');
        subplot(235)
        plot(RF_dw , 'r');
        title('Freq. Sweep')
        subplot(236)
        plot(adiabty)
        title('adiabaticity')
        pause(0.1); drawnow;
    end
    
end

seqEnd = n*modlength;

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
M = M(1:end-1,:);

%%
% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)
phi = gambar * sum(Bz * dt);



% Make nice plots here
if showMag
    
    figure(3)
    set(gcf,'Name', 'Magnetization (x,y,z)')
    set(gcf,'Position',[ 600 100 500 600 ])
    timevec = [0:(length(M)-1)]*dt;
    
    subplot(411) , plot(timevec, M(:,1))
    axis([0 timevec(end/10) -1 1]);
    subplot(412) , plot(timevec, M(:,2))
    axis([0 timevec(end/10) -1 1]);
    subplot(413) , plot(timevec, M(:,3))
    axis([0 timevec(end/10) -1 1]);
    subplot(414), plot(t,zpos);
    
    grid on
    %
    %     figure(4)
    %     plot(t,dz)
    %     title('Position of particle in time')
    %
    
end


if showSequence
    figure(1)
    set(gcf,'Name', 'Pulse sequence (x,y,z)')
    subplot(211)
    area(t, Gz)
    axis ([0 length(B1module)*dt*Npulses -5e-4 5e-4])
    subplot(212)
    plot(t, real(B1) ,'r') ,
    hold on
    plot(t, imag(B1) ,'g')
    hold off
    axis ([0 length(B1module)*dt*Npulses -5e-5 5e-5])
    
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


M = M(seqEnd,:);


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
%%
return
