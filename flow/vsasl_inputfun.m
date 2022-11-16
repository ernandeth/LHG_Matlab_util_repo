function aif = vsasl_inputfun(dt, disp, eff, R1a, tau1, tau2)
%
% function aif = vsasl_inputfun(dt, disp, eff, R1a, tau)
% 
% dt = 0.01;
% disp = 0.15;
% R1a = 1.7;
% eff = 0.9;
% tau1 = 0.1
% tau2 = 2.5;

Tmax = 3;
t = linspace(0,Tmax,Tmax/dt);

% bolus of label:
kern =  1./(exp((t-tau2)/disp^2)+1);
kern2 =  1-1./(exp((t-tau1)/disp^2)+1);
kern = kern .*kern2;

% decay of bolus:
aif = 1-2*eff*exp(-t*R1a).*kern;

%plot(t,inp)
