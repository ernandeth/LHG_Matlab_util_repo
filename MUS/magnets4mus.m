
% Ultraound calculations here:
f0 = 1e6 % target frequency of ultrasouns (Hz)
q = 1.6e-19 % single proton charge in coulombs
E = 0;  % external electric field
t = linspace(0,2e-3, 1000);
alpha = 0.1  % attenuation coefficient in water  (dB/MHz/cm)
             % in Brain it's 0.6  (why the big difference?)
c = 1480 % speed of sound in water m/s
P_peak = 5e6 % ultrasound intensity in Pascals ... might need more.
I_us = 30 % mW/cm2
I_us = 300 % W/m2

% ohm's law for sound:  P = vel * Z_acoustic
Z_acoustic = 1.51*1e6 % Rayls
v_ac = P_peak / Z_acoustic

ARF = 2*alpha*I_us / c
v_stream =  1e-6 % m/s
v_ac =  v_ac * sin(2*pi*f0*t); % m/s

% Magnetic field calculations here
ep = 8.854e-12  % farads/m
mu = 4*pi*1e-7  % Henry/m

% Toroid calculations: (T650 from micrometals)
R = 0.165/2 % Radius of toroid m
r = (0.165-0.088)/2 % cross sectional radius of toroid
a = r/2  % location of interest inside the toroid
N = round(2*pi*R / 0.02)  % Number of turns.  I want separation of 2 cm between loops.
I = 100e-3  * sin(2*pi*f0*t);  %  current through coil 

% Magnetic field in the gap of a toroid:
B = mu * N * I / (2*pi*(R+a));

% Self inductance of the toroid
L = mu * (N * r )^2 / (2*R)

%required capacitance for resonance
% 2*pi* f0 = 1 / sqrt(L* C)
C = 1 / ( L* (2*pi*f0)^2 )



% Calculate Lorentz force on particles
% F = q*(E + cross(v,B))
% if perpendicular:
F_stream = q*(E + v_stream * B)
F_ac = q*(E + v_ac .* B);
plot(F_ac)



% resulting current from movement of charged particles in fluid
% From Norton 2003 - I'm not sure about this:
% concentrations of the charges particles:
npos = 1e23
nneg = 1e23
% mobilities of the particles:
uneg = 1
upos = 1

Jresult_ac  = (npos * upos + nneg * uneg) * F_ac 
Jresult_stream  = (npos * upos + nneg * uneg) * F_stream 


