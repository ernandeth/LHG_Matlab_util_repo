function [dist1 dist2] = estimate_test

for n=1:1000
    dist1(n) = domodel(0.8);
end


for n=1:1000
    dist2(n) = domodel(0);
end
[n1, z1] = hist(dist1,50);
[n2, z2] = hist(dist2,50);

plot(z1, n1);
hold on
plot(z2, n2,'r');
hold off

cz1 = cumsum(z1)/sum(z1);
cz2 = cumsum(z2)/sum(z2);

alpha_ind = find(cz1 > 0.90);
alpha_ind = alpha_ind(1);
alpha = z1(alpha_ind);

beta_ind = find( abs(z2-alpha)<0.01)
beta = 1 - cz2(beta_ind);

return
%%

function z = domodel(effectsize)
tdim = 42;  % days
Npts = 7;
t = linspace(0,tdim, Npts);
dt = t(2) - t(1);
CNR = 50;

% self influence terms .. rate of change in (days)^-1
A = diag([0.1 0.1 0.1]);
%A = zeros(3);

% cross influence terms
B = zeros(3);
B(2,1) = effectsize;

% input influence terms

C = 0.1*eye(3);

A = A/dt; B = B/dt; C = C/dt;

% Implement a slow Input process into node 1 only
u = zeros(Npts,1);
%u(51:200) = 1 - exp(-0.03*[0:149]);
%u(201:end) = u(200)*exp(-0.01*[0:299]);

u = sin(linspace(0,pi,Npts))';
subplot(212)
plot(u)

uu = [u u u]';

%uu = randn(size(uu));


% baseline values
x = zeros(3,Npts);
x(:,1) = 0;
%x0  = [1 1 1]';  % baseline equilibrium value;
%x0  = [0 0 0]';  % baseline equilibrium value;


N = 1;
%A(:) = 0
%(:) = 0

for t=N+1:Npts
        
    %         % compute DCM for a simple model
    %     x(:,t) =  ...
    %         A * (x(:,t-1)) + ... % self-influence terms
    %         B * (x(:,t-1)) + ... % cross influence terms
    %         C * u(:, t )  ; % external input into network 1
    %
    
     x(:,t) =    (A+B)*dt * x(:,t-N) + C * uu(:, t)  ; % external input into network 1
    
end

noise = randn(size(x));
NoiseLevel = max(x(:)) / CNR;

x = NoiseLevel*noise + x;


dxdt = diff(x,1,2)/dt;



x = x';
%  x = orth(x);

X = x;
% past samples (n-1)
X1 = X(1:end-1,:);

% present samples (n)
X2 = X(2:end,:);


bhat = zeros(3,3);
RSS1 = bhat;
RSS2 = bhat;


subplot(211)
plot(X)
legend('1','2','3');
subplot(212)
drawnow

%
%% an attempt at mulitvariate regression
yy = X2;
xx = X1;
bb = pinv(xx)*yy
S = cov(yy);
varest = bb*S*bb'
tt=bb./varest



%% do EM using the model
% X2 - X1*b - cu = e ;
%
%{
cost = 10;
b = randn(3);
c = randn(3);
ctr = 0;
oldcu = zeros(size(X1));
cu  = oldcu;
u = cu;
e = cu;
ICE = 0;

while ctr<10
   %update B
    oldb = b;
    b = pinv(X1)*(X2 - u*c - e)  
    deltab = b - oldb ;
    
    %update cu  
    oldcu = cu;
    cu = X2 - X1*b - e
    %cu = cu*ICE;
    deltacu = (cu - oldcu);
    
    u = cu * pinv(c)
    plot(u); drawnow;
    
    %update C
    c = pinv(u)*(X2 - X1*b)
    
    ctr = ctr + 1;
    
    % the residuals should be more gaussian and smaller than before
    e = X2 - X1*b - u*c ;
    CE = cov(e);
    ICE = inv(CE);
    
    cost(ctr) = norm(X2 - X1*b -u*c);
    
end
plot(cost)
b
B+A

%}

% Estimate one pair at a time
for n=1:3
    for m=1:3
        %data: the present sample in this channel (nth)
        Y = X2(:,n);
        % model the previous sample in the 'influencing' channel (mth), and
        % the previous sample in this (nth) channel (to account for AR struct)
        X = [X1(:,m) X1(:,n) ones(Npts-1,1)];
        %X = [X1(:,m) X1(:,n) ];
        %X = [X1(:,m) ones(Npts-1,1)];
        
        b = pinv(X)* Y ;
        bhat(m,n) = b(1);  % influence from m to n
        
        RSS1(m,n) =  (Y - X*b)' * (Y - X*b);
        
%{      
        plot(Y,'b'), hold on,
        plot(X*b,'g' ),
        plot(X1(:,m), 'r'), hold off
        legend('original', 'fitted', 'influencing node');
        title(sprintf('Node %d -> Node %d',m, n))
        drawnow; pause
 %}       
        % Now just the nth channels's own previous history:
        X = [X1(:,n) ones(Npts-1,1)];
        b2 = pinv(X)* Y ;
        RSS2(m,n) =  (Y - X*b2)' * (Y - X*b2);


    end
end

bhat
RSS1;

t = bhat./RSS1;

z = spm_t2z(t(1,2) ,4);

% z = (RSS2 - RSS1)
fscore = (100*(RSS2-RSS1) ./ RSS2);
z = fscore(2,1);
% z =  z(2,1);
% z = tt(1,2) - tt(2,1);
% 

Bguess = 0.01*eye(3);
objfn = @(B) mymodel(B, X1, X2);
[bhat2, resnorm, residual] = lsqnonlin(objfn, Bguess, [], []);

B

return

%%
function result = mymodel(A, x1, x2data)

N = 1;
Npts = size(x1,1);
x2 = x1*A;

result = x2;

if nargin==3
    result = x2 - x2data;
end


return
