% simulation for current field imaging in magnetic ultrasound transcranial
% stimulation fields.

% calculate the magnetic field generated by a DC current through a medium

gambar = 42576000;     % gamma in Hz/T
freq = 1000;  %  Hz

duration = 2e-3; % seconds
dt = 1e-5 ; % 0.01 miliseconds per step
t = linspace(0, duration ,  duration / dt);  


current = 1e-3;  % amps

wire = [zeros(10,1) linspace(-0.05,0.05,10)', zeros(10,1)];

wire = [ 0.02*sin(linspace(0,2*pi,50))', 0.02*cos(linspace(0,2*pi,50))' zeros(50,1)];


[B Bx By Bz] = biot3d( current, wire);

subplot(321); lightbox(Bx); title('B_x (Tesla)')
subplot(323); lightbox(By); title('B_y (Tesla)')
subplot(325); lightbox(Bz); title('B_z (Tesla)'); caxis([-1e-8 1e-8])

% scale that field by a time varying waveform
wave = abs(sin(2*pi*freq* t));
subplot(322); plot(t, wave); title('current waveform')
 

% calculate the phase of the MRI signal induced by the time varying magnetic field 
phs = Bz * gambar * sum(wave)*dt;

phs = atan(tan(phs));

subplot(326); lightbox((phs)); title('MRI Phase (rads)') ; caxis([-1e-3 1e-3])
colormap(jet)

