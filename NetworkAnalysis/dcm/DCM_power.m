function DCM_power

allTPR = [];
allFPR = [];
Niter = 100;
CNR = 10;
%for CNR = [1:10]
    for r_crit = [0.5:0.1:1]
    
        [b TPR FPR ] = dosim(CNR, Niter, r_crit);
        allTPR = [allTPR TPR];
        allFPR = [allFPR FPR];
    
    end
%end

plot(allFPR, allTPR)

return




function [b TPR FPR ] = dosim(CNR, Niter, r_crit);

close all
% defaults
if nargin==0
    CNR=100;
    Niter = 100;
end

tdim = 500;
% self-influence (decay)
A = diag([-0.05 -0.05 -0.05]);
% cross influence terms
B = zeros(3);


% input influence terms
C = -diag([0.01  0.01  0.01]);

% values for the actual simulation
A1 = A; B1 =B; C1 = C;
A2 = A; B2 =B; C2 = C;

% node 1 influences node 2
B2(2,1) = -0.5;   % <--- this is the difference between the groups

% Implement a slow Input process into node 1 only
u = zeros(tdim,1);
u(51:300) = 1 - exp(-0.03*[0:249]);
u(301:end) = u(300)*exp(-0.01*[0:199]);

u = sin(linspace(0, pi,tdim))';
u = [u u u ]';

u = randn(size(u));

%u(:) = 0;

% the true function is sampled 7 times during the course of 50 days
% (~ 7 weeks)
Nsamps = 7;
samptimes = round(linspace(1, tdim, Nsamps));
t = linspace(0, 50, Nsamps);
dt = t(2)-t(1);



% compute the node activities from the model

subplot(221)
x1 = mydcm(A1, B1, C1, u, tdim)' ;
x1obs = x1(samptimes,:);
hold on;
plot(t, x1obs,'*');

subplot(222)
x2 = mydcm(A2, B2, C2, u, tdim)' ;
x2obs = x2(samptimes,:);
hold on
plot(t, x2obs,'*');

u_obs = u(samptimes);

% Niter = 50;
noiseLevel =  max(abs(x1obs(:)))/CNR;
%noiseLevel =  0;

for n=1:Niter
    noise1 = noiseLevel*randn(size(x1obs));
    
    subplot(221)
    hold on;
    plot(t, x1obs + noise1,'o');
    
    noise2 = noiseLevel*randn(size(x1obs));
    subplot(222)
    hold on;
    plot(t, x2obs + noise2, 'o');
    
    %%%%%% E-M approach
    %     [beta1, C1, u1] = dcm_em(x1obs + noise1, dt);
    %     subplot(223)
    %     x1est = mydcm( zeros(3), beta1', C1', u, tdim);
    %
    %
    %     [beta2, C2, u2] = dcm_em(x2obs + noise2, dt);
    %     subplot(224)
    %     x2est = mydcm( zeros(3), beta2', C2', u, tdim);
    %
    %     %
    
    % Now solve the model coefficients from the undersampled data
    %
    [beta1 ] = mvr_shift3(x1obs + noise1);
    [beta2 ] = mvr_shift3(x2obs + noise2);
    
    
    %beta1 = beta1(:)';
    %beta2 = beta2(:)';
    %%%%%%%%%
    
    b(n,:) = [beta1 beta2];
    %bvar(n,:) = [bvar1 bvar2];
    deltab(n,:) = ( beta1-beta2);
    
    
    
end

%t = b./bvar;
%
%return
mb = mean(b,1);
b1 = mb(1:end/2);
b2 = mb(end/2+1:end);
b1 = reshape(b1,3,3)';
b2 = reshape(b2,3,3)';

%mbvar = mean(bvar);
mbvar = var(b);
bv1 = mbvar(1:end/2);
bv2 = mbvar(end/2+1:end);
bv1 = reshape(bv1,3,3)';
bv2 = reshape(bv2,3,3)';


figure

subplot(221)
imagesc(abs(b1)); title('correlations, model 1');
caxis([0 1])
subplot(222)
imagesc(abs(b2)); title('correlations, model 2');
caxis([0 1])

subplot(223)
imagesc(abs(A1+B1)); title('True Coeffs, model 1');
caxis([0 0.1])
subplot(224)
imagesc(abs(A2+B2)); title('True Coeffs, model 2');
caxis([0 0.1])

colormap gray
% classification:
node21_1 = abs(b(:,4));   % correlation coefficients when there was no influence
node21_2 = abs(b(:,13));  % correlation coerfficients where there WAS an influence

diffy_1 = abs(b(:,4) - b(:,7));
diffy_2 = abs(b(:,13) - b(:,16));


FPR = length(find(node21_1 > r_crit))/Niter
TPR =  length(find(node21_2 > r_crit))/Niter  % <--- power

%FalsePositives = length(find(diffy_1 > 0.1))/Niter
%TruePositives =  length(find(diffy_2 > 0.1))/Niter  % <--- power


return





%%
function [Cmat] = mvr_shift3(x_obs )
% Looks for correlations between x[n-1] and x[n]
% one column at a time
% by looking at each column at a time
%



y0 = x_obs(1:end-1,:);
y1 = x_obs(2:end,:);


% for each column of the output (Y)
% estimate paramaters and their variance
% for the model   Ym = X*beta;
%
Cmat = zeros(3);
beta=[];
bvar = [];
for i=1:3
    for j=1:3
        r = corrcoef( y0(:,i), y1(:,j) );
        Cmat(i,j) = r(1,2);
        
    end
end

Cmat = Cmat(:)';

return


%%

function [B, C, u] = dcm_em(y, dt);

dydt = diff(y) / dt;
y = y(2:end,:);

dydt = dydt;
y = y;

%initial guesses:
B = -0.01 * eye(3);
C = 0.001* randn(3);
%C = diag([0.1 0.1 0.1]);

u = ones(size(y));

deltaE = 10;
oldE = zeros(size(u));
n=0;

while (deltaE > 1e-6),
    n=n+1;
    fprintf('\riteration ... %d   deltaE = %f',n, deltaE);
    
    % The model is
    %       dydt = y*B - u*C
    % calculating residuals and cov(residuals) - for diagnostics
    R = dydt - y*B - u*C;
    Sigma = cov(R);
    
    reg=0;
    
    % (1)Expectation step  :
    % expectation of uC:    E[u*C] = dydt - y*B = E[u]*C
    % so the E[u] alone:    E[u]  = (dydt - y*B)*pinv(C)
    E = (dydt - y*B ) * pinv(C)  ;
    
    % update u to its expectation
    u = E;
    
    %(2) Maximization step
    % want to find maxima of likelihood function.  For a gaussian, this happens
    % when  dydt - y*B * u*C   is minimum
    %
    % make some surrogate variables:
    alpha = dydt - y*B;   % we now minimize norm2(alpha - u*C)
    beta =  dydt - u*C;   % we now minimize norm2(beta - y*B)
    
    %reg = reg*norm(diff(u(:)));
    reg = 0.01;
    
    % maximization step:  update B and C to their least squares values
    fprintf('\ncond(y) : %f ', cond(y));
    fprintf('\ncond(u) : %f ', cond(u));
    fprintf('\ncond(B) : %f ', cond(B));
    fprintf('\ncond(C) : %f ', cond(C));
    
    B = regpinv(y, reg) * beta;
    C = regpinv(u, reg) * alpha;
    
    B = pinv(y) * beta;
    C = pinv(u) * alpha;
    
    deltaE = norm(E - oldE);
    
    oldE = E;
    
    % check this iteration:
    B
    C
    
end



return
%%






function x = mydcm(A, B, C, u, tdim)
% generate time series for three nodes
% given the connections strengths (A,B,C matrices)
% of the DCM
% the system receives a hard-coded input (u) at node 1
%

% establish some default values
if nargin==0
    % Coefficients for DCM
    % self influence terms
    A = diag([-0.05 -0.05 -0.05]);
    
    % cross influence terms
    B = zeros(3);
    B(2,1) = 0.05;
    
    % input influence terms
    C = diag([-0.01  0.01 0.01 ]);
    
    % Implement a slow Input process into node 1 only
    u = zeros(tdim,1);
    u(51:200) = 1 - exp(-0.03*[0:149]);
    u(201:end) = u(200)*exp(-0.01*[0:299]);
    u = [u u u]';
end


dt = 50/tdim;


% baseline values
x = zeros(3,tdim);
x(:,1) = 0;
%x0  = [1 1 1]';  % baseline equilibrium value;
%x0  = [0 0 0]';  % baseline equilibrium value;
N = 20;
for t=N+1:tdim
    
    
    %         % compute DCM for a simple model
    %     dx_dt =  ...
    %         A * (x(:,t-1) -x0) + ... % self-influence terms
    %         B * (x(:,t-1) -x0) + ... % cross influence terms
    %         C * u(:, t-1)  ; % external input into network 1
    %
    %     x(:,t) = x(:, t-1) + dx_dt*dt;
    
    % compute DCM for a simple model
    dx_dt =  ...
        A * mean(x(:,t-N:t-1),2 ) + ... % self-influence terms
        B * mean(x(:,t-N:t-1),2 ) + ...  % corss influence terms
        C * mean(u(:,t-N:t-1),2 ) ; %  + ...  % external input into network 1
    
    
    x(:,t) =  x(:, t-1) + dx_dt*dt;
    
end
% time units: 50 days / tdim
t = linspace(0,50,tdim);
%{
figure
subplot(211), plot(t, u,'k'); title('System input (e.g., drug)')
subplot(212),
%}
plot(t, x); title('Node Activity')
legend ('Node 1 ', 'Node 2 ', 'Node 3 ')
return


%{
function [beta, bvar] = mvr_shift2(x_obs , u_obs )
% Solves the least squares problem, on a model that is
% dy(n)/dt = B *( y(n-1) + u)
%
% by looking at each column of dy/dt at a time
%

% mean center everything:
x_obs = x_obs - repmat(mean(x_obs), [size(x_obs,1) ,1]);
u_obs = u_obs - repmat(mean(u_obs), [size(u_obs,1) ,1]);

dt = 50/length(x_obs);

y0 = x_obs(1:end-1,:);
y1 = x_obs(2:end,:);
dydt = diff(x_obs,1)/dt;


% append some ones for a DC
% dy = [ones(size(dy,1), 1) dy ];

% write the model in generic terms : Y = X*beta;
X = [ (y0+y1)/2   (u_obs(1:end-1)+u_obs(2:end))/2];
%X = [y0 ];
Y = dydt;

[t,n] = size(X);
[t,m] = size(Y);

invX = pinv(X);
sigmaY2 = var(Y);

% for each column of the output (Y)
% estimate paramaters and their variance
% for the model   Ym = X*beta;
%
beta=[];
bvar = [];
for i=1:m
    
    Ytmp = Y(:,i);
    bhat = invX * Ytmp
    sigma_bhat = diag(invX*sigmaY2(m)*eye(size(Y,1))*invX');
    
    % we only care about the second part of the model
    bhat = bhat(1:end-1);
    sigma_bhat = sigma_bhat(1:end-1);
    
    beta = [beta bhat'];
    bvar = [bvar sigma_bhat'];
end

return
%}



%{
function [beta, bvar] = mvr_shift(x_obs , u_obs )
% Solves the multivariate least squares problem, on a model that is
% delta_y(n) = B *( y(n-1) + yo + u)
%
y  = x_obs(2:end,:);
dy = diff(x_obs,1);
% dy = [ones(size(dy,1), 1) dy u_obs(2:end, :) ];
% y = [ y u_obs(2:end) ];

%{
dy_inv = pinv(dy);
beta = dy_inv*y;


C = zeros(size(beta));
bvar = zeros(size(beta));
sigma2 = var(y);

for n=1:size(C,1)
    for m=1:size(C,2)
        C = zeros(size(C));
        C(m,n) = 1;
        tmp = C'*dy_inv * sigma2(m) * dy_inv'*C;
        bvar(m,n) = tmp(m,n);
    end
end
%}
[n,p] = size(y) ;
[n,d] = size(dy);
for i = 1:n
    y_cell{i} = [kron([y(i,:)],eye(d))];
end
    
[beta, bvar, resid] = mvregress(y_cell, dy);
%beta = reshape(beta,d,d);

return
%}
