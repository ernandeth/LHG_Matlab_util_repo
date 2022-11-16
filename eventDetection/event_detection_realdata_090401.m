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
paradigm = 2   %(1 = mixed paradigm, 2 = event related paradigm)
doRebinning =0;
sigma = 0;

if paradigm == 1
	duration = 255 * TR;
	fixtime= [9 9 9 10 10 10 10 10 10 7 7 7 6 6 6 6 6 6 8 8 8];
	acttime = [10 10 10 1 1 1 1 1 1 5 5 5 1 1 1 1 1 1 10 10 10];

	D = zeros(duration / TR , 2);
	t = 1;
	for c=1:length(fixtime)
		start = t + fixtime(c)/TR;
		stop = t + fixtime(c)/TR + acttime(c)/TR;
		t=stop;
		D(start : stop-1, 2) = 1;
	end

	D(:,1) = 1;
	x = D(:,2);
	y =load('single_tdata.dat');
	y =load('sphere5mm_tdata.dat');
	y =load('BlockEvent01_onepix_tdata.dat');

	tau1 = 17;
	tau2 = 41;
end

if paradigm == 2
	duration = 306 * TR;
	fixtime= [7 9 9 9 9 9 10 10 10 10 10 10 5 5 5 5 5 5 7 7 7 7 7 7 6 6 6 6 6 6 8 8 8 8 8 8 ];
	acttime = ones(size(fixtime));

	D = zeros(duration / TR , 2);
	t = 1;
	for c=1:length(fixtime)
		start = t + fixtime(c)/TR;
		stop = t + fixtime(c)/TR + acttime(c)/TR;
		t=stop;
		D(start : stop-1, 2) = 1;
	end

	D(:,1) = 1;
	x = D(:,2);
	y =load('events_sphere_tdata.dat');
	%y =load('events_onepix_tdata.dat');
	
	tau1= 19;  % for the second data set
	tau2= 44;
end

y = mydetrend(y);
%y = y-mean(y);
y = y/max(y);
%y = y - min(y);
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



T = length(y);
p = T-10;
TP = 20;
h_err = 0;

allmTPR =[];
allmFPR =[];


H_true = HRF_mat(tau1, tau2, T);
% for s2 = 0:0.02:0.5  % iterate over noise level
% for TP = 1:20   % iterate over number of events in series
for h_err = 0;    %linspace(-2,2,5) %0.0001:0.005:0.1
% for cutoff = [0 0.1:0.1:0.4]
% for Nnbrs = [0 3:8]
% for h=[0 0.2:0.1:1]
%for sigma = 0 %[0 0.5:0.5:1]

	H = HRF_mat(tau1+h_err, tau2+h_err, T);

	if sigma ~= 0

% 		% make filtering matrix:
% 		Gkernel = normpdf( [-T:T-1], 0, sigma );
% 		Gmat =  toeplitz( [Gkernel(1); zeros(T-1,1)], Gkernel')';
% 		Gmat = Gmat(T+1:end,:);
% 
% 		% filter the data
% 		y = Gmat*y_raw;
% 		%y = y-mean(y);
% 		%y = y/max(y);
% 
% 		% also filter the model!
% 		H = HRF_mat(tau1+h_err, tau2+h_err, T);
% 		H = Gmat*H;

	end
	

	count2=1;
	for h1 = [linspace(0.0001,0.5,20) 1 3 5 ]

		fprintf('   sensitivity param:  %f  \n', h1);

		%%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
		rho = 1e-8;
		h2 = 5e-8;
		h2 = 0.0625 ; % from ROC results on 3/9/09
		%h2 = 0.1;
		%h2 = 1; % maybe?
		xhat = deconv_tv_l1(y,H,h1,h2);

		%%%%
		figure(1), subplot(211), stem(x), hold on, stem(xhat,'r'), hold off
		axis([0 250 -1.5 2])
		
		title('Event Time Course'),
		hold off, legend('True Events', 'Detected Events'); drawnow;
		figure(1), subplot(212),
		plot(y_raw), hold on, %plot(y,'r'),
		plot(H*xhat,'r'), hold off,
		title('BOLD Time Course')
		legend('Empirical Data', 'Model Fit')
		axis([0 250 -2 1])
		drawnow


		xhat(find(xhat < 1e-1))=0;
		%xhat(find(xhat < 1e-2))=0;  %threshold for calling something an activation
		xhat(find(abs(xhat)~=0))=1;

		if doRebinning
			xhatold = xhat;
			xold = x;
			
			xhat = rebin(xhat,3);
			x = rebin(x,3);
		end

		%%%%%%  count TPR and FPR  %%%%%%%%
		TPind=find(x);
		TNind=find(x==0);

		TPR(count2) = sum(x(TPind) & xhat(TPind)) / length(TPind) % sensitivity
		FPR(count2) = sum(~x(TNind) & xhat(TNind)) / length(TNind) % specificity
		count2 = count2+1;
		
		if  doRebinning
			x = xold;
			xhat = xhatold;
		end


		%		yy = H*xhat;
		%		plot(y); hold on ; plot(yy,'r');legend('data', 'estimate'); hold off
		%		drawnow

	end
	ROCarea = abs(trapz([1 FPR],[1 TPR]))
	allmTPR = [allmTPR; TPR];
	allmFPR = [allmFPR; FPR];

	
	save ROC_empirical
	% ROC curve
	figure(2)
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
	hold on ; plot(allmFPR(1,:), allmTPR(1,:), '*k');
	line([0 1], [0 1])
end

figure 
plot(allmFPR', allmTPR')
axis([0 1 0 1]);


title('ROC for mixed paradigm')
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);

%legend('Error = -2','Error = -1','Error = 0','Error = 1','Error = 2')
line([0 1], [0 1])
dofontsize(16); fatlines
