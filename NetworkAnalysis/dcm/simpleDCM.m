
% Implement a slow Input process into node 1 only
tdim = 500;
u = zeros(tdim,1);
u(51:200) = 1 - exp(-0.03*[0:149]);
u(201:end) = u(200)*exp(-0.01*[0:299]);
u = [u u u]';
dt = 1;

% Coefficients for DCM
% self influence terms
A = diag([-0.05 -0.05 -0.05]);

% cross influence terms
B = zeros(3);
B(2,1) = 0.05;

% input influence terms
C = diag([-0.01  0.01 0.01 ]);

% baseline values
x = ones(3,tdim);
x0  = [1 1 1]';  % baseline equilibrium value;
for t=2:tdim
    
  
    % compute DCM for a simple model
    dx_dt =  ...
        A * (x(:,t-1) -x0) + ... % self-influence terms
        B * (x(:,t-1) -x0) + ... % cross influence terms
        C * u(:,t-1)  ; % external input into network 1
    
    x(:,t) = x(:, t-1) + dx_dt*dt;
    
    
    
end
close all

figure
ts1 = x(1,:);
ts2 = x(2,:);
ts3 = x(3,:);
t = linspace(0,21,length(ts1));
subplot(211), plot(t, u,'k'); title('System input (e.g., drug)')
subplot(212), plot(t,ts1); title('Node Activity')
hold on
subplot(212), plot(t, ts2,'r');
subplot(212), plot(t, ts3,'g');
legend ('Node 1 Baseline', 'Node 2 Baseline', 'Node 3 Baseline')
hold off
