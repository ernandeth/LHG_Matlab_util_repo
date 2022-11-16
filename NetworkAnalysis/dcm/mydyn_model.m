% True parameters
B_true = 0.01*[1 5 3; 4 2 6; 7 8 9];
u_true = [3 2 7]';

% Simulate data
n = 50;
I = eye(3);

sig = 0.1; 
eps = sig*randn(3, n);                       % Noise N(0,sig^2)

X = rand(3, n) + eps;                               % Time 't'
for i=2:n
    X(:,i) = ( B_true) * X(:,i-1) - u_true;      % Time 't-1'
end

% Iteratively solve for B and u
u = ones(3,1);           % Initialize u (think of better initialization)
niter = 100;             % Number of iterations

for i = 1:niter
    % B update
    Y = X(:, + u * ones(1,n);
    B = I - (Y * X_new' * pinv(X_new * X_new'));
    
    % u update 
    Z = (I-B) * X_new - X_old;
    u = (Z * ones(n,1)) / n;
end



