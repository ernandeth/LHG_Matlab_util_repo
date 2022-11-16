function vel_sim_04
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
movparms = [-10, 60, 0];

% default sequence parameters:
seqparms(1) = deg2rad(10);  %flip_ang ;
seqparms(2) = 150 ; %Npulses ;
seqparms(3) = 0.5 ; %RF_duration ; % ms
seqparms(4) = 16 ;   %sweepWidth ;     % RF frequency sweep width
seqparms(5) = 0 ;   %RF_phase ;  % RF phase
seqparms(6) = 2 ;  % Total Gz_duration in ms.(includes both Gz lobes)
seqparms(7) = 1  ; %Gz_amp ; % Gauss / cm
seqparms(8) = 1;  % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(9) = 20;   % velocity target to invert

% Testing different pulse sequence parms:
seqparms(6) = 2.2 ; %Gz_duration; % duration in ms.(includes both Gz lobes)
seqparms(7) = 3  ; %Gz_amp ; % Gauss / cm
seqparms(8) = 1;  % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(1) = deg2rad(15);  %flip_ang ;
seqparms(4) = 20 ;%sweepWidth :  RF frequency sweep width (KHz)
seqparms(2) = 100 ; %Npulses ;
seqparms(3) = 0.5 ; %RF_duration ; % ms

% Testing different pulse sequence parms:  make it shorter!!  This one
% works well
seqparms(6) = 2 ; %Gz_duration; % duration in ms.(includes both Gz lobes)
seqparms(7) = 3  ; %Gz_amp ; % Gauss / cm
seqparms(8) = 1;  % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
seqparms(1) = deg2rad(30);  %flip_ang ;
seqparms(4) = 8 ;%sweepWidth :  RF frequency sweep width (KHz)
seqparms(2) = 16 ; %Npulses ;
seqparms(3) = 0.5 ; %RF_duration ; % ms
seqparms(9) = 20;   % velocity target to invert


% Going shorter to reduce off-resonance effects
% seqparms(6) = 1 ; %Gz_duration; % duration in ms.(includes both Gz lobes)
% seqparms(7) = 3  ; %Gz_amp ; % Gauss / cm
% seqparms(8) = 1;  % Gfac : fraction of refocus gradient (Gfac=1 mean balanced)
% seqparms(1) = deg2rad(40);  %flip_ang ;
% seqparms(4) = 8 ;%sweepWidth :  RF frequency sweep width (KHz)
% seqparms(2) = 16 ; %Npulses ;
% seqparms(3) = 0.5 ; %RF_duration ; % ms
% seqparms(9) = 20;   % velocity target to invert

% Trying a PCASL approach where there is no RF frequency sweep
% seqparms(1) = deg2rad(15);  %flip_ang ;
% seqparms(7) = 0.6 ;  %Gz_amp ; % Gauss / cm
% seqparms(4) = 0;      % sweep width
% seqparms(2) = 180;    % N pulses
% seqparms(8) =  1.05;  % G_fac

% sweep rate: kHz/ms = 
sweepRate = seqparms(4) / (seqparms(2) * (seqparms(6)+seqparms(3)) )



% things to try:
% flips = deg2rad([1:2:10]);
% gzamp = [2:0.4:5];
% gzdur = 4:2:16;
% npuls = 10:2:20;
% rfdur = 0.5:0.3:3;

velocities = linspace(-50,50,50);
alphas=[];
sweep_range = [2:1:64];
starting_positions = [-10:10:10];
gmax_range = [0.5:0.5:3];
flips = deg2rad(linspace(5,30,10));
init_pos = [-15:5:15];

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
% for sweep = sweep_range
%     seqparms(4) = sweep;
    


% TRy different starting positions of the spin
% for pos = init_pos
%     movparms(1) = pos;

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
        
        Mc = vel_sim(seqparms, movparms, 1, 0);
        control(n) = Mc(end,3);
    end
    
    my_inversions = [my_inversions ; inversion];    
    my_controls = [my_controls ; control];
    
    figure(18)
    hold on
    plot(velocities, inversion)
    plot(velocities, control, 'g')
  %  plot(velocities, control-inversion, 'r')
    legend('label','control','subtraction')
    title('Final Mz as function of velocity.')
    
    xlabel('Velocity')
    ylabel('M_z')
    grid on
    
    alphas = [alphas; min(inversion)];
    drawnow
end

%% Just making a nice figure to illustrate
figure(84)
subplot(313)
plot(velocities, inversion)
title('Velocity Profile of VSAI pulses.')
xlabel('Velocity (cm/s)')
ylabel('M_z')
grid on


movparms(2) = 30;

M = vel_sim(seqparms, movparms, 0, 1);


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
% 
% figure(86)
% imagesc(my_inversions - my_controls); 
% xlabel('Velocities')
% ylabel('Start position')
% ylabel('Flip Angles')
% title('Subtractions');
% colormap jet

% plot(sweep_range , alphas)
% plot(gmax_range , alphas)
% plot(rad2deg(flips) , alphas)

% movparms(2) = 50;
%
% M = vel_sim(seqparms, movparms, 0);

return



function M = vel_sim(seqparms, movparms, isControl , makeNiceFig);
%%
showSequence = 0;
showMag = 0;
slomo = 0;


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

% Manual phase correction for phase increment
MPC = 0 ;

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
    RF_phase = seqparms(5);  % RF phase
    
    % gradient pulses parameters
    Gz_duration = seqparms(6); % duration in ms.(includes all fours Gz lobes)
    Gz_amp = seqparms(7) ; % Gauss / cm
    Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm
    
    Gfac = seqparms(8);  % fraction of the gradient
    
    vel_target = seqparms(9);   % Target velocity to invert
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

% constructing the RF pulses:

% define a 0.5 ms Hanning pulse and Normalize it so the area is 1.
myHanning = hanning (RF_duration/ dt);
hlen = length(myHanning);


% Trying out other shapes of the RF envelope
% myHanning = sinc(linspace(-pi, pi,hlen))';
% myHanning(:) = 1;

% scale it to achieve flip angle - result is in Tesla
% the area required for flip_angle is  = flip / gamma
% and put it in Tesla units
B1area = flip_ang  / gam;  % area of required B1 in (Tesla * ms)
myHanning = B1area * myHanning / sum(myHanning*dt);

% Constructing the Gradient pulses:
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


% Putting RF and Gz together into a sequence:
B1seg = myHanning;

% Create an excitation module by stitching the elements together
% create some zero padding
B1pad = zeros(size(Gzseg));
Gzpad = zeros(size(myHanning));

B1module = [B1seg ; B1pad];
Gzmodule = [Gzpad ; Gzseg];
mod_length = length(B1module);
mod_duration = RF_duration + Gz_duration;

% a time vector for the duration of one module
ttmp  = linspace(0 , mod_duration, mod_length)';

% calculate the position if it moves with the target velocity and acceleration:
zpos_target = pos0 + vel_target*ttmp + 0.5*accel_target*(ttmp.^2);


% calculate the phase gained by acceleating spins during one module
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
    %     RF_phase = Mxy_phase(end) + pi;
    %RF_phase = Mxy_phase(end) - MPC;
    %RF_phase  = RF_phase + pi;
    % RF_phase = -RF_phase;
    % Gzseg = -Gzseg;
%     Gzseg(:) = 0;
%     RF_phase = 0;
    Gzseg = -Gzseg;
    %RF_phase = -RF_phase;
end

% alternative
% if isControl
%     RF_phase = -RF_phase;
% end

% design weighting function for desired acceleration profile:
% sinc for rect profile:
weights = sinc(linspace(-pi*0.6, pi*0.6,Npulses));
weights = weights * Npulses/sum(weights);
weights(:) = 1;

% frequency offset
sweepWidth = sweepWidth /dt ;  % kHz
freq_step = sweepWidth/Npulses;
RF_dw = - sweepWidth/2 ;

%RF_dw =  sweepWidth* sin(linspace(-pi/2,pi/2,Npulses));


for n = 1: Npulses
    
    FMx = cos(RF_dw  * pi * linspace(-RF_duration , RF_duration,length(myHanning))' * dt + n*RF_phase ) ;
    FMy = sin(RF_dw  * pi * linspace(-RF_duration , RF_duration,length(myHanning))' * dt + n*RF_phase ) ;
 
    RF_dw = RF_dw + freq_step ;
   
%   FMx = cos(RF_dw(n)  * pi * linspace(-RF_duration , RF_duration,length(myHanning))' * dt + n*RF_phase ) ;
%   FMy = sin(RF_dw(n)  * pi * linspace(-RF_duration , RF_duration,length(myHanning))' * dt + n*RF_phase ) ;
    
    
    
    % Hanning window envelope
    myHanning_x = myHanning .* FMx;
    myHanning_y = myHanning .* FMy;
    
    
    % RF pulses can be weighted independently
    B1seg = weights(n)*(myHanning_x  + i * myHanning_y);
    
    % update the excitation module
    B1module = [B1seg ; B1pad];
    Gzmodule = [Gzpad ; Gzseg];
    modlength = length(B1module);
    
    B1( (n-1)*modlength + 1 : n*modlength) = B1module;
    Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule;
    
    
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
    subplot(412) , plot(timevec, M(:,2))
    subplot(413) , plot(timevec, M(:,3))
    axis([0 timevec(end) -1 1]);
    subplot(414), plot(timevec,zpos);
    
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
    
    area(t, Gz)
    hold on
    plot(t, 100*real(B1) ,'r') ,
    plot(t, 100*imag(B1) ,'g')
    hold off
    %axis([0 0.1 -2 2])
    legend('Gz',  'B1x', 'B1y')
    
    
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
