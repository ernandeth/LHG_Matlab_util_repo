close all;
clear all;

count = 1
count2 = 1;
TP = 20 ; 
s2 = 0.05;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=6;
tau2=16;
h = inline('t.^(tau1-1).*exp(-t)/gamma(tau1) - 0.16*t.^(tau2-1).*exp(-t)/gamma(tau2)','tau1','tau2','t');
T = 200;
t = [0:2:2*T-2]';
hh = h(tau1,tau2,t);   hh = hh - mean(hh); hh = hh/norm(hh);
H = toeplitz( [hh(1); zeros(T-1,1)], hh')';  %  *** there was a transpose missing here ! ***
H=H(:,1:T); for i=1:T, H(:,i)=H(:,i)/norm(H(:,i)); end
H = H-ones(T,1)*mean(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count2=1;
for h1 = 0.8% linspace(0.01,3,20) 
    for count=1%:50
        count
        %%%%%% Define the event sequence %%%%%%%%%
        x=zeros(T,1);
        onsets = (T-15)* rand(TP,1) +1;
        x(round(onsets)) = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        n = sqrt(s2)*randn(T,1);
        y = H*x + n;
        y=y-mean(y);
        rho=1e-8;
        h2=5e-4;
        xhat = deconv_tv_l1(y,H,h1,h2);
        xhat(find(abs(xhat)<1e-3))=0;
        xhat(find(abs(xhat)~=0))=1;
      
        TPind=find(x==1); 
        TNind=find(x==0);

        TPR(count) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
        FPR(count) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
    end
    mTPR(count2) = mean(TPR);
    mFPR(count2) = mean(FPR);
    count2 = count2+1;
end
% ROC curve
figure; plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
line([0 1], [0 1])
title('ROC'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);