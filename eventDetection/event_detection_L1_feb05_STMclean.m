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
H = H-ones(T,1)*mean(H);

for i=1:T, H(:,i)=H(:,i)/norm(H(:,i)); end

TP = 20;
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series


for threshold =  0.0001:0.01:0.5
    for count=1:50

        %%%%%% Define the event sequence %%%%%%%%%
        x=zeros(T,1);
        % note that I'm sticking in a long activation period too.
        onsets = (T-15)* rand(TP,1) +1;
        x(round(onsets)) = 1;
        % x([10,30:50,60, 70, 90, 120, 140, 150:170])=1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % The fMRI signal y = Hx + n
        % s2=0.5;  % I'm seeing a lot of sensitivity to noise...

        n = sqrt(s2)*randn(T,1);
        y = H*x + n;


        %%%%%%%% Start Matching pursuit %%%%%%%%%%
        res = y;
        xhat = zeros(T,1);


        improvement = 1;
        old_RSS = eps;
        Nevs = 0;
        gamma=inf;
        while abs(gamma)>threshold
            e=y-H*xhat; % current residual
            C = H'*e;
            [lambda,ind] = max(abs(C));
            R=abs(H(:,ind)'*e)-threshold;
            gamma = sign(H(:,ind)'*e)*R;
            d=zeros(T,1); d(ind)=1;

            xhat = xhat + gamma*d;
            %res = y - H*xhat;

            % check to see that the detected events made a difference.
            %RSS = sum(res.^2);
            %improvement = abs(old_RSS-RSS)/old_RSS;
            %old_RSS = RSS;
            Nevs = Nevs+1;
        end

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