%function edie_spatial_100308

% Event Detection using L1,  Total variation, nonnegativity and
% majorization-minimization.
%
% version for ROC curves with spatial extent requirement
%
%close all;
clear all;
mydate = date;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;  % noise level relative to signal
% stopping rule : when we can't reduce the RSS by more than 1%

h2 = 0.01;
xdim= 5;
ydim = 5;
zdim = 5;

neighborhood = 1;
doRebinning=0;

%%%%%%%% Define HRF %%%%%%%%%%
tau1=19;
tau2=42;
T = 200;
p = T-10;
TP = 20;
h_err = 0;

allmTPR =[];
allmFPR =[];
allmTPR1 =[];
allmFPR1 =[];

allAreas=[];
allCNR = linspace(0.001,2,5);
allh1 = [ linspace(1e-4,1,20)];
%allh1 = [ linspace(1e-4,1,5)];
allh1 = [ linspace(1e-4,5,10)];
CNR = 0.5;



%H = HRF_mat(tau1+h_err, tau2, T);

count2=1;

% generate the data
%%%%%% Define the event sequence %%%%%%%%%

paradigm=1
TR = 1;


if  paradigm==1

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
	load('mixedcube.mat');
	allY=tc;
	allX = zeros(size(allY));
	for p=1:size(allY,2)
		allX(:,p) = x;
	end
	ThePix = sub2ind([xdim, ydim, zdim],3,3,3);

end
if paradigm == 2
	duration = 306 * TR;
	fixtime= [7 9 9 9 9 9 10 10 10 10 10 10 5 5 5 5 5 5 7 7 7 7 7 7 6 6 6 6 6 6 8 8 8 8 8 8 ];

	%fixtime= [0 9 9 9 9 9 10 10 10 10 10 10 5 5 5 5 5 5 7 7 7 7 7 7 6 6 6 6 6 6 8 8 8 8 8 8 ];
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

	y=load('../events_sphere_tdata.dat');

	tau1 = 18;
	tau2 = 30;

end
Npix = xdim * ydim *zdim;
T = size(y,1);
H_true = HRF_mat(tau1, tau2, T);
H = H_true;


for h1 = allh1


	y = detrend(y);
	y = y-mean(y);
	y = y/max(y);

	%%% here's the deconvolution step
	xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
	xhat_tmp = xhat;


	%%%%%%  PIXEL LEVEL: count TPR and FPR  %%%%%%%%
	threshold = 0.1;
	xhat(find((xhat)< threshold))=0;
	xhat(find((xhat)~=0))=1;



	if  0

		figure(1), subplot(211), stem(x), hold on, stem(10*xhat_tmp,'r'), hold off
		axis([0 T -1.5 2])

		title('Event Time Course'),
		hold off, legend('True Events', 'Detected Events'); drawnow;
		figure(1), subplot(212),
		plot(y), hold on, %plot(y,'r'),
		plot(H*xhat_tmp,'r'), hold off,
		title('BOLD Time Course')
		%legend('Syntehtic Data', 'Model Fit')
		axis([0 T -2 1])
		drawnow
		%pause

	end


	if doRebinning
		xhatold = xhat;
		xold = x;

		xhat = rebin(xhat,3);
		xhat(find(xhat))=1;
		
		x = rebin(x,3);
		x(find(x))=1;
	end

	% now we count for true positives, false positives
	% without the spatial requirement
	TPind=find(x==1);
	TNind=find(x==0);

	TPR = sum(x & xhat)/sum(x)
	FPR = sum(~x & xhat)/sum(~x)

	if  doRebinning
		x = xold;
		xhat = xhatold;
	end

	% with the spatial requirement


	allmTPR = [allmTPR; TPR];

	allmFPR = [allmFPR; FPR];

end

% ROC curve
figure(5)
plot( allmFPR, allmTPR, 'k'); axis([0 1 0 1]);
%line([0 1], [0 1])
title('ROC for Mixed Paradigm'); xlabel('False Positive Rate'); ylabel('True positive Rate');

dofontsize(16); fatlines
legend boxoff

axis square
drawnow


save ROC_empirical_mixed


