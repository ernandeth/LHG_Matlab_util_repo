%function edie_noise_100308

% Event Detection using L1,  Total variation and
% majorization-minimization and a negativity penalty
%
% version for ROC curves with variations in Noise level.
%
close all;
clear all;
mydate = date;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;  % noise level relative to signal
% stopping rule : when we can't reduce the RSS by more than 1%
h2 = 0.5;
NITER = 50;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=17;
tau2=42;
T = 200;
p = T-10;
TP = 20;
h_err = 0;

H_true = HRF_mat(tau1, tau2, T);
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series
% for threshold = linspace(0.1,3,20) %0.0001:0.005:0.1

allmTPR =[];
allmFPR =[];

allAreas=[];

allCNR = linspace(0.001,2,5);
%allCNR = 2;
allh1 = linspace(0.001,10,10);
allh1 = [ linspace(0.001,5,10) 10];
for CNR = allCNR  % iterate over noise level in the data
	fprintf('noise level:  %f  \n', s2);

	%H = HRF_mat(tau1+h_err, tau2, T);
	H = H_true;


	count2=1;
	for h1 = allh1

		fprintf('   sensitivity param:  %f  \n', h1);

		for count=1:NITER  % iterations for ROC samples

			%%%%%% Define the event sequence %%%%%%%%%
			x = zeros(T,1);
			onsets = (T-15)* rand(TP,1) +1;
			x(round(onsets)) = 1;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			s2=x'*H_true'*H_true*x/T/CNR;
			n = sqrt(s2)*randn(T,1);
			y = H_true*x + n;
			y = y-mean(y);
		
			%%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
			rho = 1e-8;
			xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
			xhat_tmp = xhat;

			%%%%%%  count TPR and FPR  %%%%%%%%
			threshold = 0.1*std(xhat);
			xhat(find((xhat)< threshold))=0; % 1e-3 default
			xhat(find((xhat)~=0))=1;


			TPind=find(x==1);
			TNind=find(x==0);

			TPR(count) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
			FPR(count) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
        
			if  FPR<0.2
				figure(1), subplot(211), stem(x), hold on, stem(xhat_tmp,'r--'), hold off
				axis([0 200 -1.5 2])

				title('Event Time Course'),
				hold off, legend('True Events', 'Detected Events'); drawnow;
				figure(1), subplot(212),
				plot(y), hold on, %plot(y,'r'),
				plot(H*xhat_tmp,'r'), hold off,
				title('BOLD Time Course')
				legend('Syntehtic Data', 'Model Fit')
				axis([0 200 -2 2.5])
				drawnow

				pause
			end
			
			
		end
        mTPR(count2) = mean(TPR);
        mFPR(count2) = mean(FPR);
        count2 = count2+1;
    end
    % ROC curve
    plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
    line([0 1], [0 1])
    title('ROC dependence on noise'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
    hold on
    drawnow
    allmTPR = [allmTPR; mTPR];
    allmFPR = [allmFPR; mFPR];
	
	ROCarea = abs(trapz([1 mFPR],[1 mTPR]));
	fprintf('\n   h2:  %f  ROC area: %f \n', h2, ROCarea);
	allAreas=[allAreas ROCarea];
    save ROC_noise
end

figure 
plot(allmFPR', allmTPR')
axis([0 1 0 1]); axis square


title('ROC dependence on noise'); 
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);

legend('CNR = 0.001', 'CNR = 0.5', 'CNR = 1.0','CNR =1.5', 'CNR = 2')

dofontsize(16); fatlines

