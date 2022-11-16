% Generates ASL Training Data for a machine learning based estimation
% technique. The parameters are uniformly random over a range, and Gaussian
% noise is added (mu = 0; sigma = .01) considering mtis0 = 1.
% 
%
% LHG: sep 2021
% modiifcations for VSASL labeling and new training: 2021

% number of training + testing time-series and the length of time series;
nData = 6.0*10^6;
%nData = 50;
Nframes = 150;
Y = single(zeros(nData,Nframes)); % stores the Time points;
intrv_pdf = (1.5*10^7)^-1;
sigma = 0.0; % Train data without noise - add noise in training script

% train-test split ratio
split = 0.9;

% read in the timing paramaters:
rp  = randperm(Nframes); rp=rp(1:end/2);
rp2 = randperm(Nframes); rp2=rp2(1:end/2);

timing_parms.labelcontrol           = ones(Nframes,1);
timing_parms.labelcontrol(rp)  = 0;
timing_parms.doArtSup               = ones(Nframes,1);
timing_parms.doArtSup(rp2)      = 0;

timing_parms.t_adjusts              = 0.1 *ones(Nframes,1);
timing_parms.t_delay                = 0.2 + 1.3*abs(chirp(linspace(0.75,0,Nframes), 0, 1,3.5))';
timing_parms.t_delay(2:2:end)       = timing_parms.t_delay(1:2:end);
%timing_parms.t_delay(:)       = 1;

timing_parms.t_tag                  = 0.10*ones(Nframes,1) ;
timing_parms.ArtSup_delay           = 0.15*ones(Nframes,1) ;
timing_parms.t_aq                   = 0.0329*18*ones(Nframes,1) ;
timing_parms.order                  = 1;


delays = timing_parms.t_delay;
ordr = timing_parms.order

seqDuration = sum(timing_parms.t_delay)+ ...
    sum(timing_parms.t_tag) + ...
    sum(timing_parms.t_adjusts) + ...
    sum(timing_parms.ArtSup_delay)+...
    sum(timing_parms.t_aq)
%
%figures and subtraction
doFigs = 0;
doSub = 0;
% generate Data : single run
tic

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01
parms.cbva =    0.01;
parms.kfor =    0;
parms.bat2 =    0;
parms.bat =     0.1;
parms.r1tis =   1/1.4;
parms.flip =    deg2rad(60);

parms.alpha_ti = 0.75;
parms.alpha_ai = 0.8;
parms.alpha_ts = 0.17;
    
test_signal = abs(single(  ...
            gen_signals_vs_210217(parms,...
            delays, ...
            timing_parms,...
            doFigs,doSub) ...
            + (sigma*(randn(1,Nframes))))); % LH Model one BAT
        
plot(test_signal)
drawnow
%%

f       = rand(1,nData) * 0.02;
cbv     = rand(1,nData) * 0.05;
kfor    = rand(1,nData) * 0;
bat_tis = rand(1,nData) * 0;
bat_art = rand(1,nData) * 1;
r1tis   = 0.3 + 2*rand(1,nData) ;
flip        = (pi/2)*( 0.9 + 0.2*rand(1,nData));
alpha_ti    = 0.75  *( 0.9 + 0.2*rand(1,nData));
alpha_ai    = 0.80  *( 0.9 + 0.2*rand(1,nData));
alpha_ts    = 0.17  *( 0.9 + 0.2*rand(1,nData));

theta = [
    f;
    cbv;
    kfor;
    bat_tis;
    bat_art;
    r1tis;
    flip;
    alpha_ti;
    alpha_ai;
    alpha_ts]; % LH model

theta = single(theta)';


delays = timing_parms.t_delay;

%%

        
%plot(test_signal)
%%

%%
fprintf('\nBegin generating %d time courses', nData);
%parpool('local', 'SpmdEnabled', false)
%pp =  parpool(24)
doFigs = 0;

parfor idx = 1:nData
    
    %%%% flip this for
    
    Y(idx,:) = abs(single(  ...
        gen_signals_vs_210217(genStructParmLH_VSI(theta(idx,:)),...
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

%% split into training and testing data
tic
lim = round(split*nData);
trainData = Y(1:lim,:);
trainTgts = theta(1:lim,:);

testData = Y(lim+1:end,:);
testTgts = theta(lim+1:end,:);
%%

!mkdir chirpfun/
save -v7.3 chirpfun/trainData.mat trainData
save -v7.3 chirpfun/trainTgts.mat trainTgts
save -v7.3 chirpfun/testData.mat testData
save -v7.3 chirpfun/testTgts.mat testTgts
toc

