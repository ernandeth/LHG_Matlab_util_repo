% first set of coordinates from the 1  back meta-analysis in the owen paper
%{
talcoords = [ 
    30,0,48;
    16,4,56;
    0,12,42;
    36,36,24;
    10,-58, 54;
    24,-60,52;
    10,-48,64;
    -34, -58, 42;
    42,-50,36;
    ];
ROInames = {
    'LPM1', 'LPM2', 'SMA', 'DLPFC', 'MPP1','MPP2','MPP3','IPL1', 'IPL2'};
%}

% ROI centers based on higher order N-back results in Owen paper
talcoords = [
    28 4 50;
    -26 0 52;
    -44 -2 38;
    
    -2 12 42;
    
    40 32 30;
    -36 44 20;
    -44 18 22;
    
    -30 18 6;
    32 20 4;
    
    -38 44 20;
    
    10 -66 48;
    
    -36 -50 40;
    40 -48 38;
    ];
ROInames = {
    'LPM1', 'LPM2', 'LPM3', 'SMA', 'DLPFC1','DLPFC2','DLPFC3','VLPFC1','VLPFC2','FP', 'MPP1','IPL1', 'IPL2'};


mnicoords = tal2mni(talcoords);



% now build the design matrix
Nsubjects = 24;
X = [1  1]';
X = kron(eye(Nsubjects),X);
X = [X ; X];
r = zeros(Nsubjects*4,2);
r(1:end/2, 1) = 1;
r(end/2+1:end, 2) = 1;
X = [X r];
alphareg = load('/home/data/asl/mfiles/alphaRegressor.txt');
X = [X alphareg ];

test_type = 2;

% in the first case we're testing for changes between the pre test and the
% post test... this is the same paired t-test accounting for the variance
% of the individuals, etc.  we use the same design matrix as for the whole
% group effects analysis.

if test_type==1
    contrast = zeros(1,size(X,2));
    contrast(24:25) = [-1 1];

end
    flags.header = [];
    flags.doWhiten=0;
% in the second case we test to see if, after accounting for subject's
% variance and other confounds, there is a correlation between the perfusion
% effect (amount of activation) and the performance in the task.

if test_type==2
    perform = load('/home/data/asl/mfiles/performance.dat');
    perform_reg = perform - mean(perform);
    perform_reg = perform_reg/norm(perform_reg);

    X = X(:, [1:end-3 end]);
    X = [X perform_reg];
    contrast = zeros(1,size(X,2));
    contrast(end) = 1;
end


% for c=1:size(mnicoords,1)
%     close all
%    
%     xyz = mnicoords(c,:);
%     ts = 'groupFlowEffects.img';
%     outname = ROInames{c};
%     ortho2005([],...
%         'xyz', xyz, ...
%         'tseries_file', ts, ...
%         'threshold', [],...
%         'spm_file','',...
%         'ROIsize', 12, ...
%         'ignore_origin',0,...
%         'output_name', outname, ...
%         'ROItype', 'sphere', ...
%         'interact', 1 ...
%         );
%  end
    
allROIs = [];
R = zeros(1,size(mnicoords,1));

for c=1:size(mnicoords,1)
    perform = load('/home/data/asl/mfiles/performance.dat');

    outname = ROInames{c};
    % load the activation perfusion estimates for all subjects
    f=load([outname '_rawtdata.dat']);
    
    
    % remove the guy who did terribly!
%     perform = perform( [1:11 14:61 64:end]);
%     f = f( [1:11 14:61 64:end]);
    
    
    delta_perform =  perform(end/2 +1:end) -  perform(1:end/2);
    delta_f = f(end/2+1:end) -  f(1:end/2);
    
    % remove those subjects where the flow data was outside the FOV
    % and keep track of that in the design matrix too.
    inds = find(isfinite(f));
    perform_tmp = perform(inds);
    f_tmp = f(inds);
    X2 = X(inds,:);

    
    % calculate the regression of all the random effects
    % ie - estimate the paramaters of the random effects design matrix 
    [betaCon(c) vCon(c) zscore(c)] = spmJr (f_tmp, X2, contrast,flags);
    allROIs = [allROIs f];
    load spmJr
        
    % calculate the corrected flows if we remove all the effects, except
    % fot the performance.
    f_corrected = f_tmp - X2(:,1:end-1)*beta_est(1:end-1);
    
    x = perform_tmp;
    y = f_corrected;
    figure(1)
    subplot(4,4,c); 
    plot(perform_tmp, f_corrected,'.');hold on;
    %plot(perform_tmp, f_tmp,'.r'); hold off;
    plot(perform_tmp(1:end/2), f_corrected(1:end/2),'.r'); hold off;

    delta_f = delta_f(isfinite(delta_f));
    delta_perform = delta_perform(isfinite(delta_f));
    
    % a kludge to make sure we only looo at thos e who improved
    delta_f = delta_f(delta_perform > 0.01);
    delta_perform = delta_perform(delta_perform > 0.01);
    
    Rtmp = corrcoef(delta_perform, delta_f); R(c) = Rtmp(1,2);
    
    figure(2);
    subplot(4,4,c); 
    plot(delta_perform, delta_f ,'g.');hold off;
    title(outname);
end  
R
zscore
save ROIstats betaCon vCon zscore ROInames R
save allROIs.dat allROIs -ascii
