% specify constants here:

R2 = 1/0.125;
R2a = 1/0.250;
R1 = 1/1.2;
R1a = 1/1.7;
f = 0.01; % 60 ml/min/100g
bat = 0.12;
lambda = 0.9;
dt = 0.01;

t_tag = 1.5;
t_delay = 0.2;
t_adjust = 2;
t_aq = 0.6;
flip=deg2rad(45);
nslices = 16;

selective = 0;
Tmax = t_adjust + t_tag + t_delay + t_aq;


t = linspace(0,Tmax,Tmax/dt);
% Relaxation operator:
R = diag(exp(-dt*[R2, R2, R1, R2a, R2a, R1a] ));

% perfusion/exchange operator:
% outflow:
F = diag(exp(-dt*[f/lambda , f/lambda ,f/lambda , 0,0,0]));
% inflow:
F(3,6) = 1-exp(-f*dt/lambda);

% equilibrium state:
M0 = [0 0 1 0 0 1]' ;

% initial state:
M = [0 0 1 0 0 1]';

% hi res input function:
dt2 = 1e-3;   % new temp resoultion
disp = 0.15; % dispersion
eff = 0.95;  % inversion efficiency
tau1 = bat; % leading edge
tau2 = bat+t_tag; % trailing edge
aif = vsasl_inputfun(dt2, disp, eff, R1a, tau1, tau2);
aif = aif(1:dt/dt2:end);
% use the precomputed input function:
F(6,6) = 0;
%%%

result = ones(Tmax/dt,2); 



for n = 1:(t_adjust/dt)
    M = (F * R) * M + (eye(6) - F*R) * M0;
    result(n,1) = M(3,1);
    result(n,2) = M(6,1);
end

m=1;
M(3) = 1-2*eff;  % invert the stationary spins.
for n = (t_adjust/dt):(t_adjust+t_tag+t_delay)/dt-1
    
    M = (F * R) * M + (eye(6) - F*R) * M0; 
    
    if selective==0
        if m<length(aif)
            M(6) = aif(m); m=m+1;
        end
    end
    result(n,1) = M(3,1);
    result(n,2) = M(6,1);
end

for n = (t_adjust + t_tag+t_delay)/dt : Tmax/dt
    
    M = (F * R) * M + (eye(6) - F*R) * M0; 
    M(3) = M(3) * (3)*cos(flip)^nslices;
    result(n,1) = M(3);
    result(n,2) = M(6);
end

plot(t,result)
legend('Tissue','Artery')