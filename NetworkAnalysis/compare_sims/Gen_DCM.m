close all
clc
clear all

Npts = 200;
Nnds = 4;

u = randn(Npts,Nnds);
u(u<=1.8) = 0;
u(u>1.8) = 1;

figure
for i=1:Nnds
    subplot(Nnds,2,2*i-1)
    stem(u(:,i)), title(sprintf('Inputs into Node %d',i));
end

A = -0.3*eye(Nnds);
B = zeros(Nnds);
C = eye(Nnds);

A(1,2)=0.3;
A(3,1)=0.3;


% A(2,1)=0.5;

B(1,3)=-0.2;
B(1,4)=-0.2;


% C(1,2)=1;
% C(2,3)=1;

X = zeros(Npts,Nnds);
dt = ones(1,Nnds);

for n = 2:Npts
    dXdt = ...
        A*X(n-1,:)' +  ...
        B.*(X(n-1,:)'*X(n-1,:))*ones(Nnds,1)+...
        C*u(n,:)';
    
    X(n,:) = X(n-1,:) + dXdt'.*dt;
end

th=0.2; % threshold for the event detection
for i=1:Nnds
    events(:,i)=X(:,i)>th;
end

for i=1:Nnds
    subplot(Nnds,2,2*i)
    stem(events(:,i),'g')
    hold on
    plot(X(:,i)),  title(sprintf('Activity at Node %d',i));
end

for n1=1:Nnds
    for n2=1:Nnds
        if n1==n2
            continue
        end
        atmp=events(:,n1);
        btmp=events(:,n2);
        [N(n1,n2),S(n1,n2),PNS(n1,n2)]=nec_suf(atmp,btmp);
    end
end

figure
subplot(2,3,1)
imagesc(N)
colorbar
title('Necessity')

subplot(2,3,2)
imagesc(S)
colorbar
title('Sufficiency')

subplot(2,3,4)
imagesc(A)
colorbar
title('Matrix A from DCM')

subplot(2,3,5)
imagesc(B)
colorbar
title('Matrix B from DCM')

subplot(2,3,6)
imagesc(C)
colorbar
title('Matrix C from DCM')

colormap(gray)
impixelinfo

[N;S]

