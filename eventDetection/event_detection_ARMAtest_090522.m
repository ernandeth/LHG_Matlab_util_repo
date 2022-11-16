%function event_detection_ARMAtest

% Event Detection using L1,  Total variation and
% majorization-minimization.
%
% version for ROC curves with variations in AR of the data.
%
close all;
clear all;
mydate = date;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;  % noise level relative to signal

% stopping rule : when we can't reduce the RSS by more than 1%
h1 = 0.8;

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
V = makeAR1mat(ARrho, T);

H_true = HRF_mat(tau1, tau2, T);
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series
% for threshold = linspace(0.1,3,20) %0.0001:0.005:0.1

allmTPR =[];
allmFPR =[];

allAreas=[];

for q = [0.65 0.75 0.85 0.95]% linspace(0.85,0.98,5)  % iterate over noise level in the data
	fprintf('   Autoregressive Coeff.:  %f  \n', ARrho);
    ARrho=0;
	V = makeAR1mat(ARrho,T);

	H = HRF_mat(tau1 + h_err, tau2, T);

	count2=1;

	for h1 = linspace(0.01,3,20)

		fprintf('   sensitivity param:  %f  \n', h1);

		for count=1:50  % iterations for ROC samples
			fprintf('   iteration:  %d  ', count);
			%%%%%% Define the event sequence %%%%%%%%%
			x = zeros(T,1);
			onsets = (T-15)* rand(TP,1) +1;
			x(round(onsets)) = 1;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%% generate noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %q=0.93;
            w=filter(1-q,[1 -q],sqrt(2.3/0.4*s2)*randn(T,1));
            v=sqrt(s2)*randn(T,1);
            n=w+v;
            
			%%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %n = sqrt(s2)*randn(T,1);
			y = H_true*x + n;
			y = y-mean(y);

			%%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
			rho = 1e-8;
			h2 = 5e-4;
			h2 = 0.0625 ; % from ROC results on 3/9/09
			xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
			xhat_tmp = xhat;

			
			%%%%%%  count TPR and FPR  %%%%%%%%
			xhat(find(abs(xhat)<1e-3))=0;
			xhat(find(abs(xhat)~=0))=1;

			TPind=find(x==1);
			TNind=find(x==0);

			TPR(count) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
			FPR(count) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
		
			if FPR<0.2
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

%				pause
			end
			
			
		end
		mTPR(count2) = mean(TPR);
		mFPR(count2) = mean(FPR);
		count2 = count2+1;
	end
	% ROC curve
	plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
	line([0 1], [0 1])
	title('ROC dependence on ARMA(1,1) noise'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
	hold on
	drawnow
	allmTPR = [allmTPR; mTPR];
	allmFPR = [allmFPR; mFPR];

	ROCarea = abs(trapz([1 mFPR],[1 mTPR]));
	fprintf('\n   h2:  %f  ROC area: %f \n', h2, ROCarea);
	allAreas=[allAreas ROCarea];

	%save ROC_AR1
end

figure(2)
plot(allmFPR(1:end,:)', allmTPR(1:end,:)')
axis([0 1 0 1]);


title('ROC dependence on Autoregressive Noise');
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);

legend('a = 0.65', 'a = 0.75', 'a = 0.85','a = 0.95')
line([0 1], [0 1])
%dofontsize(16); fatlines
