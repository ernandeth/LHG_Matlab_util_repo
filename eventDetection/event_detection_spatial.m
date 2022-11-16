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

allAreas=[];

%for s2 = linspace(0,0.05,5)  % iterate over noise level in the data
for s2=0.05
    fprintf('noise level:  %f  \n', s2);

    %H = HRF_mat(tau1+h_err, tau2, T);
    H = H_true;
    count2=1;
    for h1 = linspace(0.01,3,10)
    %for h1=1
%%

        fprintf('\n sensitivity param:  %f  \n', h1);
        % generate the data
        %%%%%% Define the event sequence %%%%%%%%%
        x = zeros(T,1);
        onsets = (T-15)* rand(TP,1) +1;
        x(round(onsets)) = 1;

        % track the true positives and true negatives
        TPind=find(x==1);
        TNind=find(x==0);

        %%%%%%%%%% generate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        y = H_true*x;
        y = y-mean(y);
        [allY slice] = fakeBOLDslice([xdim ydim], s2/10, y, y);
        [allX slice] = fakeBOLDslice([xdim ydim], 0, x, x);
        slice = reshape(slice,xdim,ydim);
        all_xhat = zeros(size(allY));

        for pix=1:Npix  % iterations for ROC samples
            %fprintf('\nEstimating pix n. %d \n', pix);

            y = allY(:,pix);			
%           n = sqrt(s2)*randn(T,1);
% 			y = H_true*x + n;
% 			y = y-mean(y);
            
            %%%%%%%  estimate event sequence %%%%%%%%%%%%%%%
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
        
        %all_xhat2 = zeros(size(all_xhat));
        %fprintf('\nchecking neighbors ...\n');
        %all_xhat2 = check_neighbors(all_xhat,1,xdim,ydim,1);
         
        all_xhat2 = all_xhat;

        subplot(311), imagesc(allX); title('true')
        subplot(312), imagesc(all_xhat);title('detected')
        subplot(313), imagesc(all_xhat2);title('detected and with neighbors')
        drawnow
        

        % now we count for true positives, false positives
        for pix=1:Npix
            
            x = allX(:,pix);
            xhat2 = all_xhat2(:,pix);
            
            TPind=find(x==1);
            TNind=find(x==0);

            if ~isempty(TPind)
                TPR(pix) = sum(x(TPind) & xhat2(TPind)) / length(TPind); % sensitivity
                FPR(pix) = sum(~x(TNind) & xhat2(TNind)) / length(TNind); % 1-specificity
            else
                TPR(pix)=nan;
                FPR(pix)=nan;
            end
            
            if 0%FPR<0.2            
                xhat = all_xhat(:,pix);
                sfigure(1);
                subplot(211), stem(x), hold on, stem(xhat,'r--'), stem(xhat2,'g--'), hold off
                axis([0 200 -1.5 2])

                title('Event Time Course'),
                hold off, legend('True Events', 'Detected Events'); drawnow;

                subplot(212)
                plot(allY(:,pix)), hold on, %plot(y,'r'),
                plot(H*xhat,'r'),
                plot(H*xhat2,'g'), hold off,
                title(['BOLD Time Course ' num2str(pix)])
                legend('Syntehtic Data', 'Model Fit','After spatial requisite')
                axis([0 200 -2 2.5])
                drawnow, pause(0.2)


            end


        end

        mTPR(count2) = mean(TPR(~isnan(TPR)));
        mFPR(count2) = mean(FPR(~isnan(FPR)));
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
    save ROC_spatial
end

figure
plot(allmFPR', allmTPR')
axis([0 1 0 1]);


title('ROC dependence on noise');
xlabel('False Positive Rate','Fontsize',18); ylabel('True positive Rate','Fontsize',18);

legend('CNR = Inf', 'CNR = 2.0', 'CNR = 1.0','CNR =0.67', 'CNR = 0.50')
line([0 1], [0 1])
dofontsize(16); fatlines

