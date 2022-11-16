%function event_detection_ARtest

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
NITER = 50;
TP = 10 ; % this is the number of POSTIVES
CNR = 2;% CNR=8; 0.05;  % noise level relative to signal

% stopping rule : when we can't reduce the RSS by more than 1%
h2 = 0.5;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=17;
tau2=42;
T = 200;
p = T-10;
TP = 20;
h_err = 0;
FPR = 1000;
TPR = 1000;

% initial covariance matrix for the noise (AR structure)
ARrho = 0;
%V = makeAR1mat(ARrho, T);

H_true = HRF_mat(tau1, tau2, T);

allmTPR =[];
allmFPR =[];

allAreas=[];
allh1 = [ linspace(0,5,10) 10];
allARrho =  linspace(0,0.9,5);

for ARrho = allARrho  % iterate over noise level in the data
	fprintf('   Autoregressive Coeff.:  %f  \n', ARrho);

	V = makeAR1mat(ARrho,T);

	H = HRF_mat(tau1 + h_err, tau2, T);

	count2=1;

	for h1 = allh1

		fprintf('   sensitivity param:  %f  \n', h1);

		for count=1:NITER  % iterations for ROC samples
			fprintf('   iteration:  %d  ', count);
			%%%%%% Define the event sequence %%%%%%%%%
			x = zeros(T,1);
			onsets = (T-15)* rand(TP,1) +1;
			x(round(onsets)) = 1;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			%%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			s2=x'*H_true'*H_true*x/T/CNR;
			
			[n]=ARnoise(T,ARrho,s2)*sqrt(1-ARrho^2);
            
			y = H_true*x + n;
			y = y-mean(y);
			%y = y/norm(y);

			%%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
			xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
			xhat_tmp = xhat;

			
			%%%%%%  count TPR and FPR  %%%%%%%%
			%threshold = 0.1*std(y);
			threshold = 0.1;
			xhat(find((xhat)< threshold))=0; % 1e-3 default
			xhat(find((xhat)~=0))=1;
			
			
			TPind=find(x==1);
			TNind=find(x==0);

			TPR(count) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
			FPR(count) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
		
			if  FPR<0.2
				figure(1), subplot(211), stem(x), hold on, stem(xhat_tmp,'r'), hold off
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
         		%pause
			end
			
			
		end
		mTPR(count2) = mean(TPR);
		mFPR(count2) = mean(FPR);
		count2 = count2+1;
	end
	% ROC curve
	plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
	line([0 1], [0 1])
	title('ROC dependence on autoregressive noise'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
	hold on
	drawnow
	allmTPR = [allmTPR; mTPR];
	allmFPR = [allmFPR; mFPR];

	ROCarea = abs(trapz([1 mFPR],[1 mTPR]));
	fprintf('\n   h2:  %f  ROC area: %f \n', h2, ROCarea);
	allAreas=[allAreas ROCarea];

	save ROC_AR1
end

figure
plot(allmFPR(1:end,:)', allmTPR(1:end,:)')
axis([0 1 0 1]);


title('ROC dependence on Autocorrelation');
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);

legend(num2str(allARrho(1)) , num2str(allARrho(2)), num2str(allARrho(3)),num2str(allARrho(4)),num2str(allARrho(5)),'Location','SouthEast')
legend boxoff
dofontsize(16); fatlines
axis square
