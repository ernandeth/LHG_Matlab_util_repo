% Simulation of the magnetization of particles during
% synaptic transmission (of Choline)
%
% 
clear all
close all

% SET CONSTANTS

ADC = 0.94e-9;              % in vitro choline ADC in units of m^2/s
%                             source: Nicolay et al, NMR in Biomedicine 8, p.365, 1995

%ADC = 0.94e-7;              % This one is just a test!
T = 0.05                    % experiment duration
dt = 1e-5;                  % time steps  sec
G = 1000;                   % Maximum gradient strength (G/m)
Nspins = 100;				% Number of particles to sumulate 
Nsteps = T/dt;			% Number of time steps in the simulation (add to 100msec.)
t = [1:Nsteps]*dt;			% vector of time points (sec)
T2 = 100e-3;				% Relaxation constant (in case we want to include T2 effects)
radius = 40e-9;             %compartment radius(m)

GAMMA = 4258*2*pi;	        %rad/s/G	

% The average size of the jumps in one direction is 
% based on Einstein's diiffusion equation
scale = sqrt(2*ADC*dt*Nsteps)/sqrt(Nsteps);
scale = scale/sqrt(2)      % 2 Dimensions!

% Define a gradient waveform (eg - bipolar gradients)
grad = (t >=  (Nsteps/4)*dt).*(t <= (Nsteps/2*dt)).* (-G) ... 
         + (t >= (Nsteps/2)*dt).*(t <= (3*Nsteps/4)*dt).* G ;

% compute the b value for this simulation: (units of s/m^2 )     
bval = (GAMMA*G*dt*Nsteps/4)^2 * (dt*Nsteps/4 - dt*Nsteps/12)
fprintf('\nUsing bvalue of %g s/mm^2',bval);

% Compute some expected results for restricted diffusion
so = Nspins;

D= radius^2/(2*T);
s = so*exp(-bval*D);
fprintf('\nBound max ADC: %g , min end signal: %g\n',D,s);
s = so*exp(-bval*ADC);
fprintf('\nFree ADC: %g , min end signal: %g\n',ADC,s);

% step through a number of ratios and simulate
count = 0;
last_net=[];

for ratio =1:-0.2:0 ;

    fprintf('\n Ratio: %f \n',ratio);
    str = sprintf('net%d = syn(Nsteps,Nspins,scale,grad, ratio, dt);',round(ratio*100));
    eval(str);
    str =sprintf('last=net%d(end);', round(ratio*100));
    eval(str);
    last_net = [last_net ; last];
    
end

subplot(211)
plot(t,abs(net100),'r',t,abs(net80),'b', t,abs(net40),'g')
title('Magnitude')
legend('100%','70%','40%')

subplot(212)
plot(t,angle(net100),'r',t,angle(net80),'b', t,angle(net40),'g')
title('Phase')
legend('100%','70%','40%')


figure
plot([100:-20:0],abs(last_net));
title('End Magnetization as a function of Restricted Fraction')


save sim_data.mat
fprintf('\n')
