function accel_01


% initial position, initial velocity, and acceleration 
% units are cm and seconds (for now)

in = [-5, 60, -55];
accelerations = [-200:10:200] ; 
%accelerations = -50;

inversion = zeros(size(accelerations));

for n = 1:length(accelerations)
    in(3) = accelerations(n);
    M = accel_sim(in);
    inversion(n) = M(end,3);
end

figure(18)
hold on
plot(accelerations, inversion)
title('Final Mz as function of accel.')
xlabel('Acceleration')
ylabel('M_z')

return

function M = accel_sim(in);
%%
showSequence = 0;
showMag = 0;
slomo = 0;
accel = -50 ; % in cm/s^2
vel = 55  ; % in cm/s
pos0 = -5 ;  % cm

accel_target = -120 ; % in cm/s^2

% RF pulse parameters
flip_ang = deg2rad(18);
flip_ang = deg2rad(18);
Npulses = 10;
RF_duration = 0.5 ; % ms

RF_dw = 0;     % RF frequency offset
RF_phase = 0;  % RF phase

% gradient pulses parameters
Gz_duration = 9; % duration in ms.(includes all fours Gz lobes)
Gz_amp = 3 ; % Gauss / cm
Gz_amp = Gz_amp / 1e4 ; % convert to Tesla/cm

if nargin==1
    accel = in(3) ;  
    vel = in(2) ;
    pos0 = in(1);
end
%%

% change from seconds to ms 
accel = accel * 1e-6;
vel = vel * 1e-3;
accel_target = accel_target * 1e-6;

% The usual constants
gambar = 42576;     % gamma in kHz/T
gam = gambar*2*pi;  % radians / (ms *Tesla)
mygam = gam / 1000; % radians / (s*Tesla)

T1 = 1660;  %ms
T2 = 200;   %ms

B0 = 3;    % keep things in Tesla

dt = 1e-3; % step size in milliseconds.  dt converts from samples to ms
NSTEPS = 300 / dt;  % 300 ms @  1 us steps.
t = linspace(0,300, NSTEPS)';

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

% constructing teh RF pulses:

% B1 field in mGauss
B1 = zeros(NSTEPS,1);

% define a 0.5 ms Hanning pulse and Normalize it so the area is 1.
myHanning = hanning (RF_duration/ dt);

% scale it to achieve flip angle - result is in Tesla
% the area required for flip_angle is  = flip / gamma
% and put it in Tesla units
B1area = flip_ang / gam;  % area of required B1 in (Tesla * ms)
myHanning = B1area * myHanning / sum(myHanning*dt);

% Constructing teh Gradient pulses:
% Gradient is a 1 -2 1 laplacian
Gz_length = Gz_duration / dt;

Gzseg = 1 * ones( Gz_length, 1);
Gzseg(1:end/4) = -1;
Gzseg(end*3/4+1:end) = -1;
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

ttmp  = linspace(0 , mod_duration, mod_length)';

% calculate the position of particle in one module
zpos = pos0 + vel*ttmp + accel*(ttmp.^2);

% the target velocity and acceleration do this:
zpos_target = pos0 + vel*ttmp + accel_target*(ttmp.^2);
% comment:  this screws things up although it's supposed to fix them
% 

% calculate the phase gained by acceleating spins during one module
% mygam: rad/s/Tesla
% Gz: T/cm
% dt : s

% old:  RF_phase = mygam * cumsum( zpos .* Gzmodule)*dt + pi/24;
Mxy_phase = -gam * accel_target * cumsum( ttmp.*ttmp .* Gzmodule)*dt   ;
RF_phase = Mxy_phase(end)  ; % control case:  add  +pi;

for n = 1: Npulses
    
    FMx = cos(RF_dw  * pi * linspace(-1,1,length(myHanning))' * dt + n*RF_phase ) ;
    FMy = sin(RF_dw  * pi * linspace(-1,1,length(myHanning))' * dt + n*RF_phase ) ;
    % comment - this works better with linear phase instead of quadratic.  I don't
    % understand why that is
    
    myHanning_x = myHanning .* FMx;
    myHanning_y = myHanning .* FMy;
    
    
    % an RF pulses
    B1seg = (myHanning_x  + i * myHanning_y);
    
    % update the excitation module
    B1module = [B1seg ; B1pad];
    Gzmodule = [Gzpad ; Gzseg];
    modlength = length(B1module);
    
    B1( (n-1)*modlength + 1 : n*modlength) = B1module;
    Gz( (n-1)*modlength + 1 : n*modlength) = Gzmodule;
    
end


% position of particle in cm over the whole period
dx = 0;
dy = 0;
dz = -5 + t * vel + (t.^2)*accel; % in cm

% Bz field caused by gradients expressed in Tesla.  This is a function of
% time and determines the k-space trajectory
Bz = (dx.*Gx + dy.*Gy + dz.*Gz) ;

% these are the transverse fields  due to B1 pulse (in Tesla )
Bx = real(B1);
By = imag(B1);

%%

% Put the whole thing together and run the simulation
beff = [Bx By  Bz];
M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
M = M(1:end-1,:);

%%
% calculate the amount of phase gained by the spins during the experiment.
% The units are in radians (note that the kHz cancel the ms)
phi = gambar * sum(Bz * dt);

if showMag
    
    figure(3)
    set(gcf,'Name', 'Magnetization (x,y,z)')
    timevec = [0:(length(M)-1)]*dt;
    subplot(311) , plot(timevec, M(:,1))
    subplot(312) , plot(timevec, M(:,2))
    subplot(313) , plot(timevec, M(:,3))
    
    figure(4)
    plot(t,dz)
    title('Position of particle in time')
    
    
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
drawnow

if slomo
    figure(2)
    
    for c=1:100:length(M)
        clf
        line([0 M(c,1)], [0 M(c,2)], [0 M(c,3)])
        axis([-1 1 -1 1 -1 1])
        drawnow
    end
end

return
