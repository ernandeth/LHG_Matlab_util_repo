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
h1 = 0.3;
h2 = 0.5;
xdim= 64;
ydim = 64;
Npix = xdim * ydim;
neighborhood = 1;

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
allmTPR1 =[];
allmFPR1 =[];

allAreas=[];
allCNR = linspace(0.001,2,5);
allh1 = [ linspace(1e-4,5,10)];
CNR = 1;

H = H_true;

%H = HRF_mat(tau1+h_err, tau2, T);

count2=1;

% generate the data
%%%%%% Define the event sequence %%%%%%%%%
x = zeros(T,1);
onsets = (T-15)* rand(TP,1) +1;
x(round(onsets)) = 1;

x2 = zeros(T,1);
onsets = (T-15)* rand(TP,1) +1;
x2(round(onsets)) = 1;

%%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = H_true*x;
y = y-mean(y);

y2 = H_true*x2;
y2 = y2-mean(y2);


s2=x'*H_true'*H_true*x/T/CNR;

figure(2)
[allY slice] = fakeBOLDslice([xdim ydim], s2, y, y2);
[allX slice] = fakeBOLDslice([xdim ydim], 0, x, x2);
print -dtiff fakeClusters

slice = reshape(slice,xdim,ydim);
all_xhat = zeros(size(allY));


for h1 = allh1

	for pix=1:Npix  % iterations for ROC samples
		%fprintf('\nEstimating pix n. %d \n', pix);

		y = allY(:,pix);
		
		%%% here's the deconvolution step
		xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
		xhat_tmp = xhat;


		%%%%%%  PIXEL LEVEL: count TPR and FPR  %%%%%%%%
		threshold = 0.1;
		xhat(find((xhat)< threshold))=0;
		xhat(find((xhat)~=0))=1;

		all_xhat(:,pix) = xhat;

		if  0
			x = allX(:,pix);
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

	% Now see it the positives had positive neighbors:

	all_xhat1 = zeros(size(all_xhat));
	fprintf('\nchecking neighbors ...\n');
	all_xhat1 = check_neighbors(all_xhat,1,xdim,ydim,1);


	figure(1)
	subplot(311), imagesc(allX); title('true')
	subplot(312), imagesc(all_xhat);title('detected')
	subplot(313), imagesc(all_xhat1-all_xhat);title('detected and with neighbors')
	drawnow


	% now we count for true positives, false positives
	% without the spatial requirement
	TPind=find(allX==1);
	TNind=find(allX==0);
% 
% 	TPR = sum(allX(TPind) & all_xhat(TPind)) / length(TPind) % sensitivity
% 	FPR = sum(~allX(TNind) & all_xhat(TNind)) / length(TNind) % specificity

	TPR = sum(allX(:) & all_xhat(:))/sum(allX(:))
	FPR = sum(~allX(:) & all_xhat(:))/sum(~allX(:))
	
	% with the spatial requirement

% 	TPR1 = sum(allX(TPind) & all_xhat1(TPind)) / length(TPind) % sensitivity
% 	FPR1 = sum(~allX(TNind) & all_xhat1(TNind)) / length(TNind) % specificity

	TPR1 = sum(allX(:) & all_xhat1(:))/sum(allX(:))
	FPR1 = sum(~allX(:) & all_xhat1(:))/sum(~allX(:))
	
	allmTPR = [allmTPR; TPR];
	allmTPR1 = [allmTPR1; TPR1];
	
	allmFPR = [allmFPR; FPR];
	allmFPR1 = [allmFPR1; FPR1];
	
end

% ROC curve
figure(5)
plot( allmFPR, allmTPR); axis([0 1 0 1]);
hold on
plot( allmFPR1, allmTPR1,'k'); axis([0 1 0 1]);
%line([0 1], [0 1])
title('ROC and spatial extent'); xlabel('False Positive Rate'); ylabel('True positive Rate');

dofontsize(16); fatlines
legend('N. Neighbors =0', 'N. neighbors =1', 'Location', 'SouthEast' )
legend boxoff

axis square
drawnow


    Area = abs(trapz([1 allmFPR'],[1 allmTPR']));
    Area1 = abs(trapz([1 allmFPR1'],[1 allmTPR1']));

save ROC_spatial



