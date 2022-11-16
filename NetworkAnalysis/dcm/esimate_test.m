tdim = 500;
% self influence terms
A = diag([-0.05 -0.05 -0.05]);

% cross influence terms
B = zeros(3);
B(2,1) = 0.05;

% input influence terms
C = diag([0.00  0 0]);

% Implement a slow Input process into node 1 only
u = zeros(tdim,1);
u(51:200) = 1 - exp(-0.03*[0:149]);
u(201:end) = u(200)*exp(-0.01*[0:299]);


uu = [u u u]';
dt = 1;


% baseline values
x = zeros(3,tdim);
x(:,1) = 1;
%x0  = [1 1 1]';  % baseline equilibrium value;
%x0  = [0 0 0]';  % baseline equilibrium value;
for t=2:tdim
    
    
    %         % compute DCM for a simple model
    %     dx_dt =  ...
    %         A * (x(:,t-1) -x0) + ... % self-influence terms
    %         B * (x(:,t-1) -x0) + ... % cross influence terms
    %         C * u(:, t-1)  ; % external input into network 1
    %
    %     x(:,t) = x(:, t-1) + dx_dt*dt;
    
    % compute DCM for a simple model
    dx_dt =  ...
        A * (x(:,t-1) ) + ... % self-influence terms
        B * (x(:,t-1) ); %  + ... % cross influence terms
    
    x(:,t) =  x(:, t-1) + dx_dt*dt + C * uu(:, t) ; % external input into network 1
    
    
end

figure

plot(x')

dxdt = diff(x,1,2)/dt;

X = [x' u];
X = X(1:end-1,:);
Y = dxdt(1,:)';
bhat = pinv(X)* Y


