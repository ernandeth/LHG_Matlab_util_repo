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
p=190;
t = [0:2:2*T-2]';
hh = h(tau1,tau2,t);   hh = hh - mean(hh); hh = hh/norm(hh);
H = toeplitz( [hh(1); zeros(T-1,1)], hh')';  %  *** there was a transpose missing here ! ***
H=H(:,1:T-10); 
for i=1:p, H(:,i)=H(:,i)/norm(H(:,i)); end
H = H-ones(T,1)*mean(H);
TP = 20;
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series


for threshold = linspace(0.1,3,20) %0.0001:0.005:0.1
    for count=1:50

        %%%%%% Define the event sequence %%%%%%%%%
        x=zeros(p,1);
        % note that I'm sticking in a long activation period too.
        onsets = (p-15)* rand(TP,1) +1;
        x(round(onsets)) = 1;
        %x([60:70])=1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % The fMRI signal y = Hx + n
        % s2=0.5;  % I'm seeing a lot of sensitivity to noise...
        
        n = sqrt(s2)*randn(T,1);
        y = H*x + n;
        y=y-mean(y);

        %%%%%%%% Start Matching pursuit %%%%%%%%%%
        y0 = y;
        x0 = y(1:p);% zeros(T,1); % Step 1


        improvement = 1;
        old_RSS = eps;
        Nevs = 0;
        rho=1e-8;
        h2=5e-4;
        while improvement > threshold/50000     % We need to find a way to pick this number
            delta_xt=diff(x0); 
            w1=h2/2./sqrt(delta_xt.^2+rho);
            W1=diag(w1);
            D=diag(ones(p,1))+diag(-ones(p-1,1),1); 
            D=D(1:p-1,:);
            w=threshold/2./sqrt(x0.^2+rho);
            W=diag(w);
            A=H'*H+W+D'*W1*D;
            x1=A\(H'*y);
           %figure(1); plot(x1); pause
            res = y - H*x1;
            % check to see that the detected events made a difference.
            RSS = sum(res.^2);
            improvement = abs(old_RSS-RSS)/old_RSS;
            old_RSS = RSS;
            Nevs = Nevs+1;
            x0=x1;
        end
        xhat=x1;
        xhat(find(abs(xhat)<1e-3))=0;
        xhat(find(abs(xhat)~=0))=1;
        allNevs(count) = Nevs;

        %Display HRF %
        %figure(1); set(gca,'Fontsize',12);
        %plot(t,hh); xlabel('t [s]','Fontsize',20); ylabel('AU','Fontsize',20); title('HRF','Fontsize',20);
        %Display the fMRI signal
        %figure(3); set(gca,'Fontsize',12);
        %
        %     plot(t,y,'k'); xlabel('t [s]','Fontsize',20); ylabel('AU','Fontsize',20); title('fMRI signal','Fontsize',20);
        %     hold on
        %     %why doesn't the signal return to baseline???
        %     stem(t,x); xlabel('t [s]','Fontsize',12); ylabel('AU','Fontsize',20); title('Estimated event sequence','Fontsize',20);
        %     stem(t,xhat, '--r'); xlabel('t [s]','Fontsize',12); ylabel('AU','Fontsize',20); title('Estimated event sequence','Fontsize',20);
        %     drawnow
        %     hold off
        
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