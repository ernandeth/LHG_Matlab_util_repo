clear
%%%%%%%%%%% Set Figure Font Sizes %%%%%%%%%%%%%%%%%%%
set(0, 'DefaultAxesLineWidth', 6.0)
set(0, 'DefaultTextFontSize', 16)
set(0, 'DefaultTextFontWeight', 'normal')
set(0, 'DefaultAxesFontSize', 16)
%set(0, 'DefaultAxesFontWeight', 'bold')
set(0, 'DefaultAxesFontWeight', 'normal')
set(0, 'DefaultLineMarkerSize', 8)
set(0, 'DefaultLineLinewidth', 4.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ran = 0;
rng(ran)
global schedflag
global nobat_tis
%%%%%%%%%%% SET OPTIMIZATION SPECS %%%%%%
Nframes = 800;
Nparms = 7; % 5 for two, 4 for single compartment
Ntrials = 25;

pr = 5 ; % number of coeffs

max_c = 2; % maximum coeff value
min_c = 0; % minimum coeff value
intv = 0.2;

W = ones(Nparms,Nparms);
W(:,4) = 2;
W(4,:) = 2;

scan_time = 540;

nobat_tis = 1;

one_bat = boolean([1 0 1 1 1 1 1]); % no bat_tis in this case
if nobat_tis
    W = W(one_bat',one_bat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wt = diag(W);
%%%% for preloaded label-control
% load('../t_parms_opt/240s/minPoly_revrand_TR.mat')
% isLabel = timing_parms_min.isLabel;
% clear timing_parms_min

%%%% for generated label-control
iL = (rand(1,Nframes));iL(iL<0.4) = 0;iL(iL>=0.8) = -1;iL(iL<0.8 & iL>=0.4) = 1;
isLabel = iL;



% parm order :: cbva bat_tis bat_art f eta r1tis flip
f=         60 /6000;
cbva =      0.01 ;
eta=   0.5 ;
bat_tis = 1.2;
bat_art = 0.6;
r1tis =     1/1.4  ;
flip =      90*pi/180 ;


params = [(0.2*cbva) + ((1.4*cbva) - (0.2*cbva)) * rand(1,Ntrials);...
        (0.3*bat_tis) + ((1.4*bat_tis) - (0.3*bat_tis)) * rand(1,Ntrials);...
    (0.3*bat_art) + ((1.4*bat_art) - (0.3*bat_art)) * rand(1,Ntrials);...
    (0.2*f) + ((1.4*f) - (0.2*f)) * rand(1,Ntrials);...
    (0.5*eta) + ((1.4*eta) - (0.5*eta)) * rand(1,Ntrials);...
    (0.3*r1tis) + ((3*r1tis) - (0.3*r1tis)) * rand(1,Ntrials);...
    (0.6*flip) + ((1.3*flip) - (0.6*flip)) * rand(1,Ntrials)];

% params = reshape(params,Nparms,1,Ntrials); 
% norm_val = sqrt(params).*sqrt(permute(params,[2 1 3]));
% norm_val = 1./norm_val;

params = reshape(params,Nparms,1,Ntrials); 
norm_val = bsxfun(@times, sqrt(params), sqrt(permute(params,[2 1 3])));
norm_val = 1./norm_val;

if nobat_tis
    norm_val = norm_val(one_bat',one_bat,:);
end


% W =diag(sqrt(wt)); % for two compartment (without T1 and Flip) %[1.5 1.5 0.25 1.5 0.25 ]
% W =diag(sqrt([1 1 1 1 ])); % for single compartment
del = 0.025; % for 'h' in finite differences



% len = 100000; % number of random schemes tested
% ObjFun = zeros(len,length(pr)); % array of objective function values
% tims = zeros(len,length(pr)); % array of total times;
crlb_mins = zeros(Nparms,length(pr));

obj_min = Inf*ones(1,length(pr)); % initializer for minimum objective function
ptr = 1; % for the position of the random distro being tested
doSub = 0;% 1 (Yes) or 0 (No)
sigLen = ( (doSub*Nframes/2) + ((1-doSub)*Nframes) );

pr_ptr = 0; % pointer to the current number of interpolation points
range = min_c:intv:max_c;
rangelen= length(range);
dims = length(range)*ones(1,pr);

% constant times per acquisition
t_delay_i = 0.0550;
t_adjust_i = 0.053;
t_aq_i = 0.0329;
tag_min = 0;

t_delay = t_delay_i*ones(1,Nframes);
t_adjust = t_adjust_i*ones(1,Nframes);
t_aq = t_aq_i*ones(1,Nframes);


timing_parms_fix.t_aq = t_aq(:);
timing_parms_fix.t_delay = t_delay(:);
% timing_parms_fix.t_tag = t_tag(comb,:)';
timing_parms_fix.t_adjust = t_adjust(:);
timing_parms_fix.Nlabel_group = 1;
timing_parms_fix.order = 1; % use the default order (usually control-tag)
timing_parms_fix.isLabel = isLabel;



tau = t_delay_i + t_adjust_i + t_aq_i + tag_min; % constant itmes per acquisition
tic

%%
for npts_pr = pr
    
    pr_ptr = pr_ptr+1;
    
    
    idx_len = rangelen^npts_pr;
    
     
    coeffs = cell(idx_len,npts_pr);
    obj  = zeros(1,idx_len); % preallocate array for cost functions 
    
    parfor comb = 1:idx_len
%%
        % using schedule:
            % Scheme(3) Polynomial Optimization 
            
            % preallocating the del_s(theta)
                D{comb} = zeros(Nparms,sigLen);
           for trials = 1:Ntrials     
            % calculate del_s(theta)0
            for par = 1:Nparms
                D{comb}(par,:)=differ_opt_rand_v2_7(params(:,:,trials),del,par,...
                    scaleTags(min_c,npts_pr,intv,rangelen,comb,Nframes,tau,scan_time,tag_min),... % Just the scaled tagging schedule
                    timing_parms_fix,doSub); %'differ_opt' for opt test
            end
                if nobat_tis
                    D{comb} = D{comb}(one_bat',:);
                end
                obj(comb) = obj(comb) + (sum(sum(W.*(sqrt(abs(inv(D{comb}*D{comb}')/10^4)).*norm_val(:,:,trials)).*W))/Ntrials); % Cost function evaluation
                D{comb} = [];
           end
    end
    
    
    [obj_min(pr_ptr),ix] = min(obj);
    timing_parms_min(pr_ptr) = timing_parms_fix;
    timing_parms_min(pr_ptr).t_tag = scaleTags(min_c,npts_pr,intv,rangelen,ix,Nframes,tau,scan_time,tag_min);
    c_mins = (min_c*ones(1,npts_pr))+ (intv*(my_ind2sub(rangelen*ones(1,npts_pr),ix)-1));
    
    params2 = params;
    if nobat_tis
        params2 = params(one_bat',:,:);
    end
    
    fprintf('Optimization complete for for %d loops! \n',npts_pr)
    
end
%%
toc;
CRLB_opt = zeros(Nparms,length(pr));
if nobat_tis
    CRLB_opt = zeros(Nparms-1,length(pr));
end
for trials = 1:Ntrials
    for i = 1:length(pr)
        CRLB_opt(:,i) = CRLB_opt(:,i) + (100*(optCRLB_rand7(params(:,:,trials),timing_parms_min(i),Nframes,doSub)./params2(:,:,trials))/Ntrials);
    end
end
CRLB_opt

%%
fname = ['7pnb_min' schedflag 'rand' '_len_' num2str((pr)) '_' num2str(tag_min) 'max' '_' num2str(max_c) '_' num2str(scan_time) '_rewtLH2_' num2str(wt') '_' num2str(Nframes) '.mat'];
save(fname,'timing_parms_min','obj_min','CRLB_opt','c_mins');

mailprefs;
sendmail('anishl@umich.edu','Experiment Completed',['The results are stored in file: ',fname]);