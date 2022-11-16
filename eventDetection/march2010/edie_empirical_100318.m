%function edie_spatial_100308

% Event Detection using L1,  Total variation, nonnegativity and
% majorization-minimization.
%
% version for ROC curves with spatial extent requirement
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
h2 = 0;
xdim= 5;
ydim = 5;
zdim = 5;

doRebinning=1;
neighborhood = 1;

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
allh1 = [ linspace(1e-4,2,20)];
allh1 = [ linspace(1e-4,1,10)];
CNR = 0.5;



%H = HRF_mat(tau1+h_err, tau2, T);

count2=1;

% generate the data
%%%%%% Define the event sequence %%%%%%%%%

paradigm=3
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
	load('events_cube.mat');
	%tc=load('../events_sphere.dat');

	zdim=4;

	allY=tc;
	allX = zeros(size(allY));
	for p=1:size(allY,2)
		allX(:,p) = x;
	end

	ThePix = sub2ind([xdim, ydim, zdim],3,3,2);
	%ThePix = sub2ind([xdim, ydim, zdim],4,2,3);
	tau1 = 18;
	tau2 = 30;

end

if paradigm == 3
	xdim= 3;
	ydim = 3;
	zdim = 3;
	
	duration = 610 * TR;
	fixtime= 30*ones(1,20);
	acttime = 0.5*ones(size(fixtime));

	D = zeros(duration / TR , 2);
	t = 1;

	for c=1:length(fixtime)
		start = t + fixtime(c)/TR;
		stop = t + fixtime(c)/TR + acttime(c)/TR;
		t=stop;
		D(floor(start) : ceil(stop-1), 2) = 1;
	end

	D(:,1) = 1;
	x = D(:,2);
	load('030314ms_vis_timecourse/cube_030314ms.mat');
	% tc=load('../events_sphere.dat');

	allY=tc;
	allX = zeros(size(allY));
	for p=1:size(allY,2)
		allX(:,p) = x;
	end

	ThePix = sub2ind([xdim, ydim, zdim],2,2,2);

	tau1 = 20;
	tau2 = 30;

end


Npix = xdim * ydim *zdim;
T = size(allY,1);
H_true = HRF_mat(tau1, tau2, T);
H = H_true;

xold = x;

if doRebinning

	x = rebin(x,3);
	x(find(x))=1;

	allX = zeros(size(allY));
	allX = allX(1:3:end-1,:);
	for p=1:size(allY,2)
		allX(:,p) = x;
	end

end

% filter stuff:
sigma=1;
		
for h1 = allh1

	for pix=1:Npix  % iterations for ROC samples
		%fprintf('\nEstimating pix n. %d \n', pix);

		y = allY(:,pix);

		Gkernel = normpdf([-20:20], 0, sigma);
		y = y-mean(y);
		tmp = conv(Gkernel, y );
		y = tmp(21:end-20);

		y = detrend(y);
		
		y = y/max(y);

		%%% here's the deconvolution step
		xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
		xhat_tmp = xhat;



		%%%%%%  PIXEL LEVEL: count TPR and FPR  %%%%%%%%
		threshold = 0.1;
		threshold = 0.05;
		%threshold = 0.1*std(y);
		xhat(find((xhat)< threshold))=0;
		xhat(find((xhat)~=0))=1;

	if  1

			figure(1), subplot(211), stem(xold), hold on, stem(0.5*xhat,'r'), hold off
			axis([0 T -1.5 2])

			title('Event Time Course'),
			hold off, legend('True Events', 'Detected Events'); drawnow;
			figure(1), subplot(212),
			plot(y), hold on, %plot(y,'r'),
			plot(H*xhat_tmp,'r'), hold off,
			title('BOLD Time Course')
			legend('Syntehtic Data', 'Model Fit')
			axis([0 T -2 1])
			drawnow
			%pause
	end

	if doRebinning
		xhatold = xhat;

		xhat = rebin(xhat,3);
		xhat(find(xhat))=1;
	end


		
	all_xhat(:,pix) = xhat;
	% Now see it the positives had positive neighbors:
	end
	all_xhat1 = zeros(size(all_xhat));
	fprintf('\nchecking neighbors ...\n');
	all_xhat1 = check_neighbors(all_xhat,1,xdim,ydim,zdim);

	figure(1)
	subplot(311), imagesc(allX); title('true')
	subplot(312), imagesc(all_xhat);title('detected')
	subplot(313), imagesc(all_xhat1-all_xhat);title('eliminated positives (blue)')
	drawnow


	% now we count for true positives, false positives
	% without the spatial requirement
	TPind=find(allX==1);
	TNind=find(allX==0);

	TPR = sum(allX(:, ThePix) & all_xhat(:, ThePix))/sum(allX(:, ThePix));
	FPR = sum(~allX(:, ThePix) & all_xhat(:, ThePix))/sum(~allX(:, ThePix));

	% with the spatial requirement

	TPR1 = sum(allX(:, ThePix) & all_xhat1(:, ThePix))/sum(allX(:, ThePix));
	FPR1 = sum(~allX(:, ThePix) & all_xhat1(:, ThePix))/sum(~allX(:, ThePix));

	allmTPR = [allmTPR; TPR];
	allmTPR1 = [allmTPR1; TPR1];

	allmFPR = [allmFPR; FPR];
	allmFPR1 = [allmFPR1; FPR1];

	[allmFPR allmTPR]
	[allmFPR1 allmTPR1]

end

% ROC curve
figure(5)
plot( allmFPR, allmTPR); axis([0 1 0 1]);
hold on
plot( allmFPR1, allmTPR1,'k'); axis([0 1 0 1]);
%line([0 1], [0 1])
title('ROC for E.R. Paradigm'); xlabel('False Positive Rate'); ylabel('True positive Rate');

dofontsize(16); fatlines
legend('N. Neighbors =0', 'N. neighbors =1', 'Location', 'SouthEast' )
legend boxoff

axis square
drawnow


save ROC_empirical_rebin


