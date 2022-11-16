% using the algorithm on empirical data.  time coursese from a sphere

% Event Detection using L1,  Total variation, nonnegativity and
% majorization-minimization.
%
% version for ROC curves with spatial extent requirement
%

clear all;
mydate = date;

count = 1;
count2 = 1;
TP = 10 ; % this is the number of POSTIVES
s2 = 0.05;  % noise level relative to signal
% stopping rule : when we can't reduce the RSS by more than 1%
h2 = 0.1;
xdim= 5;
ydim = 5;
zdim = 5;

doRebinning=0;
neighborhood = 1;

%%%%%%%% Define HRF %%%%%%%%%%

T = 200;
p = T-10;
TP = 20;
h_err = 0;

allmTPR =[];
allmFPR =[];
allmTPR1 =[];
allmFPR1 =[];

allAreas=[];
allh1 = [ linspace(1e-4,3,20)];



%H = HRF_mat(tau1+h_err, tau2, T);

count2=1;
pix=1;
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
    tc = load('../BlockEvent01_sphere_tdata.dat');
    allY=tc;
    allX = zeros(size(allY));
    for p=1:size(allY,2)
        allX(:,p) = x;
    end
    ThePix = sub2ind([xdim, ydim, zdim],3,3,3);

    tau1=17;
    tau2=41;
    tau1=16;
    tau2=40;
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
    tc = load('../events_sphere_tdata.dat');

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
if paradigm==3
    duration = 310 * TR;
    fixtime= [...
        9 9 9 9 9 9 ...
        10 10 10 10 10 10  ...
        5 5 5 5 5 5 ...
        7 7 7 7 7 7 ...
        6 6 6 6 6 6 ...
        8 8 8 8 8 8 ];
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
    tc = load('sphere12_tdata.dat');

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
y = tc;

Npix = xdim * ydim *zdim;
T = size(allY,1);
H_true = HRF_mat(tau1, tau2, T);
H = H_true;
pix=1;
ThePix=1;
xold = x;

% filter stuff:
sigma=2;

for h1 = allh1


    % 		Gkernel = normpdf([-20:20], 0, sigma);
    % 		y = y-mean(y);
    % 		tmp = conv(Gkernel, y );
    % 		y = tmp(21:end-20);

    y = detrend(y);
    y = y/max(y);
    y = y-mean(y);


    %%% here's the deconvolution step
    xhat = deconv_tv_l1_nonneg(y,H,h1,h2);
    xhat_tmp = xhat;

    %%%%%%  PIXEL LEVEL: count TPR and FPR  %%%%%%%%
    threshold = 0.1;
    threshold = 0.05;
    threshold = 0.1*std(xhat);
    xhat(find((xhat)< threshold))=0;
    xhat(find((xhat)~=0))=1;


    if doRebinning
        xhatold = xhat;

        xhat1 = rebin(xhat,1);
        xhat1(find(xhat1))=1;

        x1 = rebin(x,1);
        x1(find(x1))=1;
    end


    % now we count for true positives, false positives
    % without the spatial requirement
    TPind=find(x==1);
    TNind=find(x==0);

    TPR = sum(x & xhat)/sum(x);
    FPR = sum(~x & xhat)/sum(~x);

    allmFPR = [allmFPR; FPR];
    allmTPR = [allmTPR; TPR];

    [allmFPR allmTPR]

    if doRebinning
        % with the rebinning
        TPind1=find(x1==1);
        TNind1=find(x1==0);

        TPR1 = sum(x1 & xhat1)/sum(x1);
        FPR1 = sum(~x1 & xhat1)/sum(~x1);

        allmTPR1 = [allmTPR1; TPR1];
        allmFPR1 = [allmFPR1; FPR1];


        [allmFPR1 allmTPR1]
    end

    if  FPR<0.2

        figure(1), subplot(211),
        stem(xold), hold on, stem(0.5*xhat,'r'), hold off
        %stem(x1), hold on, stem(0.5*xhat1,'r'), hold off
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


end

% ROC curve
figure(5)
plot( allmFPR, allmTPR); axis([0 1 0 1]);

if doRebinning
    hold on
    plot( allmFPR1, allmTPR1,'k'); axis([0 1 0 1]);
    legend('No re-binning','Re-binning=1','Location', 'SouthEast' )
end

%line([0 1], [0 1])
title('ROC for E.R. Paradigm');
if paradigm==1
    title('ROC for Mixed Paradigm');
end
xlabel('False Positive Rate'); ylabel('True positive Rate');

dofontsize(16); fatlines
legend boxoff

axis square
drawnow

Area = abs(trapz([1 allmFPR'],[1 allmTPR']));
Area1 = abs(trapz([1 allmFPR1'],[1 allmTPR1']));

save ROC_empirical_rebin


