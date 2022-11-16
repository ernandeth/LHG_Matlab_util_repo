%function event_detection_Noisetest

% Event Detection using L1,  Total variation and
% majorization-minimization.
%
% version for ROC curves with variations in Noise level.
%
addpath ('~hernan/matlab/eventDetect')

close all;
clear all;
mydate = date;
randn('seed',1);
count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
% noise level relative to signal
% stopping rule : when we can't reduce the RSS by more than 1%
h1 = 0.8;
%%%%%%%% Define HRF %%%%%%%%%%
tau1=17;
tau2=42;
T = 200;
p = T-10;
TP = 20;
h_err = 0;

H_true = HRF_mat(tau1, tau2, T);
figure(1); set(gca,'Fontsize',12);
plot(H_true(:,1)); xlabel('Time [s]','Fontsize',18); ylabel('AU','Fontsize',18);
%figure(1); print -depsc C:\event_detection\latex\fig1.eps

x = zeros(T,1);
onsets = (T-15)* rand(TP,1) +1;
x(round(onsets)) = 1;
CNR=0.25;
CNR = 2;
s2=x'*H_true'*H_true*x/T/CNR;
n = sqrt(s2)*randn(T,1);
y = H_true*x + n;
figure(2); plot(x); set(gca,'Fontsize',12);
xlabel('Time [s]','Fontsize',18); ylabel('x_t','Fontsize',18); axis([0 200 0 2]);
%figure(2); print -depsc C:\event_detection\latex\fig2.eps


figure(3); plot(y); set(gca,'Fontsize',12);
xlabel('Time [s]','Fontsize',18); ylabel('y_t','Fontsize',18);
SNR=x'*H_true'*H_true*x/s2/T
logSNR=10*log10(SNR)
%figure(3); print -depsc C:\event_detection\latex\fig3.eps


if(1)
    % for s2 = 0:0.02:0.5  % iterate over noise level
    % for TP = 1:20   % iterate over number of events in series
    % for threshold = linspace(0.1,3,20) %0.0001:0.005:0.1

    allmTPR =[];
    allmFPR =[];

    allAreas=[];

    for h2=linspace(1,0,5);  % iterate over TV penalty
        fprintf('h2 parm:  %f  \n', h2);

        %H = HRF_mat(tau1+h_err, tau2, T);
        H = H_true;
        count2=1;
        %for h1 = linspace(0,6,20)  % iterate over L1 penalty
        for h1 = [0 linspace(0.001,1,10) 2:10] % iterate over L1 penalty

            fprintf('   h1 param:  %f  \n', h1);

            for count= 1:100  % iterations for ROC samples

                %%%%%% Define the event sequence %%%%%%%%%
                x = zeros(T,1);
                onsets = (T-15)* rand(TP,1) +1;
                x(round(onsets)) = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                n = sqrt(s2)*randn(T,1);
                y = H_true*x + n;
                y = y-mean(y);

                %%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
                rho = 1e-8;
                %h2 = 5e-4;
                %h2 = 0.0625 ; % from ROC results on 3/9/09
                xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
                xhat_tmp = xhat;


                %%%%%%  count TPR and FPR  %%%%%%%%
                %threshold = 0.1*std(xhat);
				threshold = 0.1;
                xhat(find((xhat)< threshold))=0; % 1e-3 default
                xhat(find((xhat)~=0))=1;

                TPind=find(x==1);
                TNind=find(x==0);

                TPR(count) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
                FPR(count) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
                
            end
            
            if  0%  FPR<0.2
                figure(1), subplot(211), stem(x), hold on,
                stem(xhat_tmp,'r--'), stem(xhat,'g--'), hold off
                axis([0 200 -1.5 2])

                title('Event Time Course'),
                hold off, legend('True Events', 'Detected Events'); drawnow;
                figure(1), subplot(212),
                plot(y), hold on, %plot(y,'r'),
                plot(H*xhat_tmp,'r'), hold off,
                title('BOLD Time Course')
                legend('Syntehtic Data', 'Model Fit')
                axis([0 200 -2 1])
                drawnow


            end
            mTPR(count2) = mean(TPR);
            mFPR(count2) = mean(FPR);
            fprintf('\t TPR = %f \t FPR = %f \n', mTPR(count2), mFPR(count2));
            count2 = count2+1;
        end
        % ROC curve
        plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
        line([0 1], [0 1])
        title('ROC dependence on noise'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
        hold on
        drawnow
        allmTPR = [allmTPR; mTPR];  % each new h2 is a row, h1 are columns
        allmFPR = [allmFPR; mFPR];



        ROCarea = abs(trapz([1 mFPR],[1 mTPR]));
        fprintf('\n   h2:  %f  ROC area: %f \n', h2, ROCarea);
        allAreas=[allAreas ROCarea];
        save ROC_testh1h2
    end

    allh1 = [0 linspace(0.001,1,10) 2:10];
    allh2 = linspace(1,0,5);
    subplot(221) , plot(allh1,allmTPR); title('TPR vs. h1')
    subplot(222) , plot(allh1,allmFPR);title('FPR vs. h1')
    subplot(223) , plot(allh2,allmTPR);title('TPR vs. h2')
    subplot(224) , plot(allh2,allmFPR);title('FPR vs. h2')

    figure(4)
    plot(allmFPR', allmTPR')
    axis([0 1 0 1]);
	axis square

 
	text(0.6, 0.5, 'h_2 value:', 'Fontsize', 16)
	legend(num2str(allh2(1)),num2str(allh2(2)),num2str(allh2(3)),...
		num2str(allh2(4)),num2str(allh2(5)),...
		'Location', 'SouthEast') 
	legend boxoff
	
    title('ROC dependence on h_1 for different h_2');
    xlabel('False Positive Rate'); ylabel('True positive Rate');

    %legend('CNR = Inf', 'CNR = 2.0', 'CNR = 1.0','CNR = 0.7', 'CNR = 0.5')
    %line([0 1], [0 1])
    dofontsize(16); fatlines
end
