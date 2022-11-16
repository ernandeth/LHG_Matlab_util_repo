function X=gen_net

Npts = 1000;

X = zeros(Npts,5);
% u = zeros(Npts,5);

% define an external influence on one of the nodes:
% inds = floor(linspace(1,Npts,50));
% u(inds,1) = 1;
% 
% inds2 = floor(linspace(4,Npts-Npts/5,3));
% for n=1:100
%     u(inds2 + n, 4) = 1;
% end

u = randn(Npts,5);
u(u<2.3) = 0;
u(u>0) = 1;

% subplot(5,1,1)
% stem(u(:,1))
% subplot(5,1,2)
% stem(u(:,2))
% subplot(5,1,3)
% stem(u(:,3))
% subplot(5,1,4)
% stem(u(:,4))
% subplot(5,1,5)
% stem(u(:,5))

ind=find(u(:,1)==1)
% u(ind(5)-5:ind(5)+5,1)=1;
u(ind(5)-5:ind(5)+5,2)=1;


% figure,
% subplot(5,1,1)
% stem(u(:,1))
% subplot(5,1,2)
% stem(u(:,2))
% subplot(5,1,3)
% stem(u(:,3))
% subplot(5,1,4)
% stem(u(:,4))
% subplot(5,1,5)
% stem(u(:,5))


%u(:,2) = 0;

% save u u

dt = ones(1,5);

% self-coefficient (decay)
A = -0.1 * ones(1,5);

% xternal influence coefficient
B =  ones(1,5);
B(4)=0;             % No external stimuli into 4.
                    % This will make it so that node 3 is necessary for node 4

% Direct influence coefficients - Sufficiency
C = zeros(5,5);
% for n=1:4
%     C(n,n+1) = 0.5;
% end
% C(2,3) = 0;
C(1,2) = 0.3;  % Node 1 is sufficient for node 2
C(3,4) = 0.3;  % node 3 is sufficient for node 4


% calculate the  changes in activity at each time step
for n=2:Npts
    dXdt = ...
        X(n-1,:) .* A  +  ...
        u(n,:) .*B  +  ...
        X(n-1,:)* C;
    
    X(n,:) = X(n-1,:) + dXdt .* dt ;
end

X = X + 0.1*randn(size(X));