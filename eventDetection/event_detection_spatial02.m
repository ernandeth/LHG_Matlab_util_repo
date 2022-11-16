%function event_detection_spatial

% Event Detection using L1,  Total variation and
% majorization-minimization.
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
h1 = 0.8;
xdim= 16;
ydim = 16;
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

for s2 = linspace(0,0.1,5)  % iterate over noise level in the data
%for s2=0.05
    fprintf('noise level:  %f  \n', s2);

    %H = HRF_mat(tau1+h_err, tau2, T);
    H = H_true;
    count2=1;

    %for h1=1
    %%

    fprintf('\n sensitivity param:  %f  \n', h1);
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
    
	figure(2)
    [allY slice] = fakeBOLDslice([xdim ydim], s2, y, y2);
    [allX slice] = fakeBOLDslice([xdim ydim], 0, x, x2);
    slice = reshape(slice,xdim,ydim);
    all_xhat = zeros(size(allY));
    
    for h1 = linspace(0.01,3,10)
        for pix=1:Npix  % iterations for ROC samples
            %fprintf('\nEstimating pix n. %d \n', pix);

            y = allY(:,pix);

            rho = 1e-8;
            h2 = 5e-4;
            h2 = 0.0625 ; % from ROC results on 3/9/09
            %h2=0; % removing TV penalty

            %%% here's the deconvolution step
            xhat = deconv_tv_l1(y,H,h1,h2);
            xhat_tmp = xhat;


            %%%%%%  pix TPR and FPR  %%%%%%%%
            xhat(find(abs(xhat)<1e-3))=0;
            xhat(find(abs(xhat)~=0))=1;

            all_xhat(:,pix) = xhat;

        end

        % Now see it the positives had positive neighbors:

        all_xhat2 = zeros(size(all_xhat));
        fprintf('\nchecking neighbors ...\n');
        all_xhat2 = check_neighbors(all_xhat,1,xdim,ydim,1);


        figure(1)
        subplot(311), imagesc(allX); title('true')
        subplot(312), imagesc(all_xhat);title('detected')
        subplot(313), imagesc(all_xhat2);title('detected and with neighbors')
        drawnow


        % now we count for true positives, false positives
        % without the spatial requirement
        TPR = all_xhat & allX;
        TPR = sum(TPR(:)) / sum(allX(:));

        FPR = all_xhat & ~allX;
        FPR = sum(FPR(:)) / sum(~allX(:));

        % with the spatial requirement
        TPR1 = all_xhat2 & allX;
        TPR1 = sum(TPR1(:)) / sum(allX(:));

        FPR1 = all_xhat2 & ~allX;
        FPR1 = sum(FPR1(:)) / sum(~allX(:));


        mTPR(count2) = TPR;
        mFPR(count2) = FPR;
        mTPR1(count2) = TPR1;
        mFPR1(count2) = FPR1;

        count2 = count2+1;
    end
    allmTPR = [allmTPR; mTPR];
    allmTPR1 = [allmTPR1; mTPR1];
    allmFPR = [allmFPR; mFPR];
    allmFPR1 = [allmFPR1; mFPR];

    % ROC curve
    figure(5)
    plot( mFPR, mTPR,'o', mFPR,mTPR); axis([0 1 0 1]);
    hold on
    plot( mFPR1, mTPR1,'*k', mFPR1,mTPR1,'k'); axis([0 1 0 1]);
    line([0 1], [0 1])
    title('ROC dependence on noise'); xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);
    hold on
    drawnow

    ROCarea = abs(trapz([1 mFPR],[1 mTPR]));
    fprintf('\n   h2:  %f  ROC area: %f \n', h2, ROCarea);
    allAreas=[allAreas ROCarea];
    save ROC_spatial
end
CNR = 1./sqrt(linspace(0,0.1,5))


