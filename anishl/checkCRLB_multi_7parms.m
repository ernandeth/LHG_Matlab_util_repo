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
ran = 10;
rng(ran)

global nobat_tis

%%%%%%%%%%% SET OPTIMIZATION SPECS %%%%%%
Nframes = length(timing_parms_min.t_tag);
Nparms = 7; % 5 for two, 4 for single compartment
Ntrials = 50;

pr = 5 ; % number of coeffs

max_c = 2; % maximum coeff value
min_c = 0; % minimum coeff value
intv = 0.2;

scan_time = 600;

W = ones(Nparms,Nparms);
W(:,4) = 2;
W(4,:) = 2;

nobat_tis = boolean(1);

one_bat = boolean([1 0 1 1 1 1 1]); % no bat_tis in this case


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% parm order :: cbva bat_tis bat_art f eta r1tis flip
f=         60 /6000;
cbva =      0.01 ;
eta=   0.02 ;
bat_tis = 2.0;
bat_art = 2.0;
r1tis =     1/1.4  ;
flip =      90*pi/180 ;


params = [(0.2*cbva) + ((1.4*cbva) - (0.2*cbva)) * rand(1,Ntrials);...
        (0.3*bat_tis) + ((1.4*bat_tis) - (0.3*bat_tis)) * rand(1,Ntrials);...
    (0.3*bat_art) + ((1.4*bat_art) - (0.3*bat_art)) * rand(1,Ntrials);...
    (0.2*f) + ((1.4*f) - (0.2*f)) * rand(1,Ntrials);...
    (0.0*eta) + ((1.0*eta) - (0.0*eta)) * rand(1,Ntrials);...
    (0.3*r1tis) + ((3*r1tis) - (0.3*r1tis)) * rand(1,Ntrials);...
    (0.6*flip) + ((1.3*flip) - (0.6*flip)) * rand(1,Ntrials)];

params = reshape(params,Nparms,1,Ntrials); 
norm_val = bsxfun(@times, sqrt(params), sqrt(permute(params,[2 1 3])));
norm_val = 1./norm_val;

if Ntrials == 1
    params = [cbva bat_tis bat_art f eta r1tis flip]';
end

params2 = params;
if nobat_tis
    norm_val = norm_val(one_bat',one_bat,:);
    params2 = params(one_bat',:,:);
end

% W =diag(sqrt(wt)); % for two compartment (without T1 and Flip) %[1.5 1.5 0.25 1.5 0.25 ]
del = 0.025; % for 'h' in finite differences

doSub = 0;% 1 (Yes) or 0 (No)
%%

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