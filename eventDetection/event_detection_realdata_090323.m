% Read real data time courses and extract events from them
% Event Detection using L1,  Total variation and
% majorization-minimization.
%
% version for ROC curves with variations in Noise level.
%
close all;
clear all;

TR = 1;
% notes - we collected 300 seonds of data but the paradigm ended after 255.
% clipped off the last 45 ffrom the end of the rarun_01.nii file

duration = 255 * TR;
fixtime= [9 9 9 10 10 10 10 10 10 7 7 7 6 6 6 6 6 6 8 8 8];
acttime = [10 10 10 1 1 1 1 1 1 5 5 5 1 1 1 1 1 1 10 10 10];

D = zeros(duration / TR , 2);
t = 0;
for c=1:length(fixtime)
	start = t + fixtime(c)/TR
	stop = t + fixtime(c)/TR + acttime(c)/TR
	t=stop
	D(start : stop, 2) = 1;
end

D(:,1) = 1;
x = D(:,2);
y =load('single_tdata.dat');
%y =load('sphere5mm_tdata.dat');

y = y-mean(y);
y = y/max(y);
y_raw = y;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;  % noise level relative to signal
% stopping rule : when we can't reduce the RSS by more than 1%
h1 = 0.8;

% Non local filter defaults:
Nnbrs = 5;
h = 1;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=6 * 2;
tau2=16 *2;
T = length(y);
p = T-10;
TP = 20;
h_err = 0;

allmTPR =[];
allmFPR =[];


H_true = HRF_mat(tau1, tau2, T);
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series
% for h_err = 0;linspace(3,0,5) %0.0001:0.005:0.1
% for cutoff = [0 0.1:0.1:0.4]
% for Nnbrs = [0 3:8]
% for h=[0 0.2:0.1:1]
for sigma = [0 1:4]
	
	if sigma ~= 0
		%y = smoothdata(y_raw,TR,cutoff,8);
		
		%y = nlfilt_090205(y_raw, Nnbrs, h);
		
		Gkernel = normpdf([-20:20], 0, sigma);
		tmp = conv(Gkernel, y_raw-mean(y_raw) );
		y = tmp(21:end-20);
		
		y = y-mean(y);
		y = y/max(y);
        
      
	end

	H = HRF_mat(tau1+h_err, tau2+h_err, T);
	count2=1;
	for h1 = linspace(0.01,3,10)

		fprintf('   sensitivity param:  %f  \n', h1);

		%%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
		rho = 1e-8;
		h2 = 5e-4;
		h2 = 0.0625 ; % from ROC results on 3/9/09
        h2 = 0.2;
		xhat = deconv_tv_l1(y,H,h1,h2);


		%%%%%%  count TPR and FPR  %%%%%%%%
		xhat(find(abs(xhat)<1e-2))=0;
		xhat(find(abs(xhat)~=0))=1;

		TPind=find(x==1);
		TNind=find(x==0);

		TPR(count2) = sum(x(TPind) & xhat(TPind)) / length(TPind); % sensitivity
		FPR(count2) = sum(~x(TNind) & xhat(TNind)) / length(TNind); % specificity
		count2 = count2+1;

		%		yy = H*xhat;
		%		plot(y); hold on ; plot(yy,'r');legend('data', 'estimate'); hold off
		%		drawnow

	end
	allmTPR = [allmTPR; TPR];
	allmFPR = [allmFPR; FPR];

	save ROC_onepix_smoothfiltered
	% ROC curve

	plot( FPR, TPR,'.', FPR,TPR); axis([0 1 0 1]);
	line([0 1], [0 1])
	title('ROC in real data'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
	hold on
	drawnow

end

if 0
	% make figs
	figure

	load ROC_onepix_smoothfiltered
	subplot(221), plot(allmFPR', allmTPR'), title('Parks McLellan Filter')
	hold on ; plot(allmFPR(1,:), allmTPR(1,:), 'k');
	line([0 1], [0 1])

	load ROC_onepix_NLfiltered
	subplot(222), plot(allmFPR', allmTPR'), title('Non-Local Filter')
	hold on ; plot(allmFPR(1,:), allmTPR(1,:), 'k');
	line([0 1], [0 1])

	load ROC_onepix_Gaussfiltered
	subplot(223), plot(allmFPR', allmTPR'), title('Gaussian Filter')
	hold on ; plot(allmFPR(1,:), allmTPR(1,:), 'k');
	line([0 1], [0 1])
end
