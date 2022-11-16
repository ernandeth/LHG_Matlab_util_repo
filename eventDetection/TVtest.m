close all;
clear all;

count = 1
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;
% stopping rule : when we can't reduce the RSS by more than 1%
% threshold = 0.01;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=6;
tau2=16;

h = inline('t.^(tau1-1).*exp(-t)/gamma(tau1) - 0.16*t.^(tau2-1).*exp(-t)/gamma(tau2)','tau1','tau2','t');
T = 200;
t = [0:2:2*T-2]';
hh = h(tau1,tau2,t);   hh = hh - mean(hh); hh = hh/norm(hh);
H = toeplitz( [hh(1); zeros(T-1,1)], hh')';  %  *** there was a transpose missing here ! ***
%H = H-ones(T,1)*mean(H);
%H=eye(T);
%%%%%% Define the event sequence %%%%%%%%%
x=zeros(T,1);
x(100:120)=1;
n = sqrt(s2)*randn(T,1);
y = H*x + n;
y=y-mean(y);
x0 = y;
Nevs = 0;
figure(1); plot(y); pause
threshold=2;
rho=0.00001;
for t=1:50     % We need to find a way to pick this number    
    delta_xt=diff(x0); 
    w=threshold/2./sqrt(delta_xt.^2+rho);
    W=diag(w);
    D=diag(ones(T,1))+diag(-ones(T-1,1),1); D=D(1:T-1,:);
    A=H'*H+D'*W*D;
    x1=A\y;
    figure(1); plot(x1); pause
    res = y - H*x1;
    % check to see that the detected events made a difference.
   
    x0=x1;
end
xhat=x1;
        
        