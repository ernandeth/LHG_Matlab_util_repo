function edie_HRFerr

% Event Detection using L1,  Total variation and
% majorization-minimization.
%
% version for ROC curves with variations in HRF.
%
close all;
clear all;
mydate = date;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
h2 = 0.5;
NITER=50;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=17;
tau2=42;
T = 200;
p = T-10;
TP = 20;

CNR = 2;

H_true = HRF_mat(tau1, tau2, T);
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series
% for threshold = linspace(0.1,3,20) %0.0001:0.005:0.1

allmTPR =[];
allmFPR =[];
allh_err = linspace(-2,2,5);
allh1 = linspace(0.001,10,10);
allh1 = [ linspace(0.001,5,10) 10];

for h_err = allh_err  % iterate over errors in the HRF model
	fprintf('    HRF error:  %f  \n', h_err);

	H = HRF_mat(tau1+h_err, tau2, T);

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

			if 0% FPR<0.2
				figure(1), subplot(211), stem(x), hold on, stem(xhat_tmp,'r--'), hold off
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

		end
		mTPR(count2) = mean(TPR);
		mFPR(count2) = mean(FPR);
		count2 = count2+1;
	end
	% ROC curve
	plot( mFPR, mTPR,'.', mFPR,mTPR); axis([0 1 0 1]);
    axis square
    %line([0 1], [0 1])
	title('ROC dependence on errors in HRF'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
	hold on
	drawnow
	allmTPR = [allmTPR; mTPR];
	allmFPR = [allmFPR; mFPR];
    
    
    save ROC_HRF
end

for count=1:size(allmTPR,1)
    ROCarea(count) = abs(trapz([1; allmFPR(:, count);0],[1; allmTPR(:,count);0]));
end

figure
plot(allmFPR', allmTPR')
axis([0 1 0 1]);


title('ROC dependence on Errors in the HRF');
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate')

legend('Error = -2','Error = -1','Error = 0','Error = 1','Error = 2',...
	 'Fontsize',18, 'Position','SouthEast');
% line([0 1], [0 1])
dofontsize(16); fatlines
legend boxoff


