% Generates ASL Training Data for a machine learning based estimation
% technique. The parameters are uniformly random over a range, and Gaussian
% noise is added (mu = 0; sigma = .01) considering mtis0 = 1.
% 
%Model: Su et al.
%
%- Anish Lahiri, 1 Aug 2017
%
% LHG: sep 2021
% modiifcations for VSASL labeling and new training: 2021

sigma = 0.0/sqrt(1); % noise std deviation

% train-test split ratio
split = 0.9;

% read in the timing paramaters:

formatSpec = '%f';
fileID = fopen('t_delays.txt','r');
timing_parms.t_delay = fscanf(fileID,formatSpec);
fileID = fopen('t_adjusts.txt','r');
timing_parms.t_adjusts = fscanf(fileID,formatSpec);
fileID = fopen('doArtSuppression.txt','r');
timing_parms.doArtSup =  fscanf(fileID,formatSpec);
fileID = fopen('isVelocitySelective.txt','r');
timing_parms.labelcontrol =  fscanf(fileID,formatSpec);

% some preset constant values:
Nframes = length(timing_parms.t_delay);
timing_parms.t_tag  = 0.1*ones(Nframes,1) ;
timing_parms.ArtSup_delay  = 0.15*ones(Nframes,1) ;
timing_parms.t_aq  = 0.0329*18*ones(Nframes,1) ;
timing_parms.order = 1;

%{
% scheme = 'minLinInterp600LH2_700mod4_20'; % 'perlin' ,'gaussian', 'sinusoid','uniform','minPoly_rand','minPoly_revrand','minLinInterp360','minPoly360LH500_20','minPoly360LH2_500_20'
scheme = 'optVSIMRF_vsson_201asup2';
% scheme = 'subopt600LH2_700_20';
filePath = ['../t_parms_opt/600_VSI/' scheme '.mat']
load(filePath)

iL = timing_parms.isLabel;
tmp = iL;
tmp(iL == 1) = 0;
tmp(iL == 0) = 1;
timing_parms.isLabel = tmp;
scheme = ['flpd_' scheme];
timing_parms = timing_parms; % optimal timing parameters for 'scheme'

%}
%%
%Distribution:
%mean1,mean2,sigma1,sigma2
dispflg = 'dd';
vsn = 12; dp.vsn = vsn;
disp =Inf; if (disp == Inf); dispflg = ''; end 
dp.mean_f_1 = 0.001; dp.mean_f_2 = 0.014; dp.sigma_f_1 = 20; dp.sigma_f_2 = 20;
dp.mean_cbv_1 = 0.001; dp.mean_cbv_2 = 0.014; dp.sigma_cbv_1 = 20; dp.sigma_cbv_2 = 20;
dp.mean_bat_1 = 0.3; dp.mean_bat_2 = 0.6; dp.sigma_bat_1 = 20; dp.sigma_bat_2 = 20;
dispflg
% if timing_parms.order == 2
%     timing_parms.order = 1
% end

ordr = timing_parms.order
%%
%figures and subtraction
doFigs = 0;
doSub = 0;

% fixed parms
parms.mtis0 =     1 ;
parms.r1blood = 1/1.7;
% parms.flip =      90*pi/180 ;

% mean values for parameter ranges and extremities
f=     60 /6000 ;   f_min = 0.1*f ;             f_max = 0.015 ;
cbva =     0.01 ;   cbva_min = 0.00*cbva ;       cbva_max= 0.015 ;
bat =   1.2 ;   bat_tis_min = 0.1 ; bat_tis_max = 0.6;
r1tis =   1/1.4 ;   r1tis_min = (1/3) ;  r1tis_max = (1/.3);
flip = 90*pi/180;   flip_min = 0.8*flip;        flip_max = 1.2*flip; % 0.4, 1.4 2019-03-07 AL

% % vsn 7
% alpha_ti = 0.66;    alpha_ti_min = 0.9*alpha_ti; alpha_ti_max = 1.1*alpha_ti;
% alpha_ai = 0.8057;    alpha_ai_min = 0.9*alpha_ai; alpha_ai_max = 1.1*alpha_ai;
% alpha_ts = 0.17;    alpha_ts_min = 0.9*alpha_ts; alpha_ts_max = 1.1*alpha_ts;% 0.77


% % vsn 8
% % updates from Luis 3/13/20
% alpha_ti = 0.75;    alpha_ti_min = 0.9*alpha_ti; alpha_ti_max = 1.1*alpha_ti;
% alpha_ai = 0.9;    alpha_ai_min = 0.9*alpha_ai; alpha_ai_max = 1.1*alpha_ai;
% alpha_ts = 0.17;    alpha_ts_min = 0.9*alpha_ts; alpha_ts_max = 1.1*alpha_ts;% 0.77

% vsn 11
% updates from Luis 3/13/20
alpha_ti = 0.75;    alpha_ti_min = 0.9*alpha_ti; alpha_ti_max = 1.1*alpha_ti;
alpha_ai = 0.8057;    alpha_ai_min = 0.9*alpha_ai; alpha_ai_max = 1.1*alpha_ai;
alpha_ts = 0.17;    alpha_ts_min = 0.9*alpha_ts; alpha_ts_max = 1.1*alpha_ts;% 0.77


% number of training + testing time-series and the length of time series;
nData = 6.0*10^6;
intrv_pdf = (1.5*10^7)^-1;
Y = single(zeros(nData,Nframes)); % stores the Time points;


%%

f_pdf = genpdf_2gausunif(dp.mean_f_1, dp.mean_f_2, dp.sigma_f_1, dp.sigma_f_2, f_min:f_max*intrv_pdf:f_max);
cbva_pdf = genpdf_2gausunif(dp.mean_cbv_1,dp.mean_cbv_2,dp.sigma_cbv_1,dp.sigma_cbv_2,cbva_min:cbva_max*intrv_pdf:cbva_max);
bat_pdf = genpdf_2gausunif(dp.mean_bat_1,dp.mean_bat_2,dp.sigma_bat_1,dp.sigma_bat_2,bat_art_min:bat_art_max*intrv_pdf:bat_art_max);

% f_s = f_min + (f_max-f_min)*rand(1,nData);
f_s = randpdf(f_pdf,f_min:f_max*intrv_pdf:f_max,[1 nData]);
% f_s = randpdf(genpdf_2gausunif(0.0005,0.014,0.003,0.003,f_min:f_max*intrv_pdf:f_max),f_min:f_max*intrv_pdf:f_max,[1 nData]);

% cbva_s = cbva_min + (cbva_max-cbva_min)*rand(1,nData);
cbva_s = randpdf(cbva_pdf,cbva_min:cbva_max*intrv_pdf:cbva_max,[1 nData]);
% cbva_s = randpdf(genpdf_2gausunif(0.002,0.014,0.003,0.003,cbva_min:cbva_max*intrv_pdf:cbva_max),cbva_min:cbva_max*intrv_pdf:cbva_max,[1 nData]);

%%%% flip this for
kfor_s = kfor_min + (kfor_max-kfor_min)*rand(1,nData); % LH model
% eta_s = eta_min + (eta_max-eta_min)*rand(1,nData); % PS Model
%%%%
bat_tis_s = bat_tis_min + (bat_tis_max-bat_tis_min)*rand(1,nData);
bat_art_s = randpdf(bat_pdf,bat_art_min:bat_art_max*intrv_pdf:bat_art_max,[1 nData]);
r1tis_s = r1tis_min + (r1tis_max-r1tis_min)*rand(1,nData);
flip_s = flip_min + (flip_max -flip_min)*rand(1,nData);

alpha_ti_s = alpha_ti_min + (alpha_ti_max -alpha_ti_min)*rand(1,nData);
alpha_ai_s = alpha_ai_min + (alpha_ai_max -alpha_ai_min)*rand(1,nData);
alpha_ts_s = alpha_ts_min + (alpha_ts_max -alpha_ts_min)*rand(1,nData);

%%%% flip this for
parm_s = [
    f_s;
    cbva_s;
    kfor_s;
    bat_tis_s;
    bat_art_s;
    r1tis_s;
    flip_s;
    alpha_ti_s;
    alpha_ai_s;
    alpha_ts_s]; % LH model
% parm_s = [f_s;cbva_s;eta_s;bat_tis_s;bat_art_s;r1tis_s;flip_s]; % PS model
%%%%


theta = single(parm_s');


delays = timing_parms.t_delay;
j = sqrt(-1);
%%

% FOR 234 SCHEDULE ONLY! take out otherwise
% timing_parms.doArtSup(1:4)=-1;

% try to figure out discrepancy between model and data
% timing_parms.doArtSup(10:end)= -1;

% generate Data
tic
mdl = '';

test_signal = abs(single(  ...
            gen_signals_vs_210217(genStructParmLH_VSI(parm_s(:,1000)),...
            delays,...
            timing_parms,...
            doFigs,doSub) ...
            + (sigma*(randn(1,Nframes))))); % LH Model one BAT
        
plot(test_signal)
%%
        
fprintf('\nBegin generating %d time courses', nData);
%parpool('local', 'SpmdEnabled', false)
%pp =  parpool(24)
parfor idx = 1:nData
    
    %%%% flip this for
    
    Y(idx,:) = abs(single(  ...
        gen_signals_vs_210217(genStructParmLH_VSI(parm_s(:,idx)),...
        delays,...
        timing_parms,...
        doFigs,doSub) ...
        + (sigma*(randn(1,Nframes))))); % LH Model one BAT
    
    if mod(idx,1e3)==0
        fprintf('\ngenerating time course number ... %d', idx);
    end
end

fprintf('\n ....Finished generating %d time courses', nData);

toc
%%
%Y = Y./repmat(Y(:,1),1,Nframes);
tic
lim = round(split*nData);
%% split into training and testing data
trainData = Y(1:lim,:);
trainTgts = theta(1:lim,:);

testData = Y(lim+1:end,:);
testTgts = theta(lim+1:end,:);
%%
toc


save -v7.3 trainData.mat trainData
save -v7.3 trainTgts.mat trainTgts
save -v7.3 testData.mat testData
save -v7.3 testTgts.mat testTgts

% Not sure what this does ... ?

% n_r = 91; %91
% n_c = 109; %109
% sigma2 = 0.003;

% % testData = testData(1:n_r*n_c,:);
% % testData = testData + sigma*(randn(size(testData)));
% % testTgts = testTgts(1:n_r*n_c,:);
% % filename2 = ['/export/anishl/nNetPredData/Data_randist_'  num2str(scheme) '_realsigma' num2str(sigma) '_ordr' num2str(ordr) 'vsn' num2str(vsn) dispflg '.mat'];
% % mask = ones(n_r,n_c);
% % save(filename2,'testData','testTgts','mask','dp','-v7.3');

toc
%}
