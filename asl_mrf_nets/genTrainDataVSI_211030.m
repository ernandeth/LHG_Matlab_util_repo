function genTrainData_VSS_211030

% Generates ASL Training Data for a machine learning based estimation
% technique. The paramT2lossrs are uniformly random over a range, and Gaussian
% noise is added (mu = 0; sigma = .01) considering mtis0 = 1.
%
%
% LHG: sep 2021
% modiifcations for VSASL labeling and new training: 2021

% number of training + testing time-series and the length of time series;
nData = 6.0*10^6;
%nData = 50;

Nframes = 401;
%figures and subtraction
doFigs = 0;
doSub = 0;

Y = single(zeros(nData,Nframes)); % stores the Time points;
intrv_pdf = (1.5*10^7)^-1;
sigma = 0.0; % Train data without noise - add noise in training script

% train-test split ratio
split = 0.9;
best=0;
% Make timing paramaters with differnt randomizations. Look for best RMS
% difference in signal change due to CBF change
for t=1:1000
    
    % random distribution of control and label pulses
    rp  = randperm(Nframes);
    rp=rp(1:end/2);
    % random distribution of arterial suppression (1/4 of the time)
    rp2 = randperm(Nframes);
    rp2=rp2(1:end/4);
    % random distribution of delays:
    rp3 = randperm(Nframes);
    %
    timing_parms.label_type         = 'FTVSI-sinc';
    %timing_parms.label_type         = 'BIR8inv';
    timing_parms.labelcontrol       = ones(Nframes,1);
    timing_parms.labelcontrol(rp)   = 0;
    timing_parms.doArtSup           = ones(Nframes,1);
    timing_parms.doArtSup(rp2)      = 0;
    
    
    % new schedule : 1/11/22
    timing_parms.t_delay = 0.3 + 0.1*(chirp(linspace(1,0,Nframes), 0, 1, 13))';
    timing_parms.t_adjusts = 0.4 +  0.2*(chirp(linspace( 0,1, Nframes), 0, 1,13))';
    
    %{
    % Anish schedule:
    tmp = load('timing_files_201/t_adjusts.txt');
    timing_parms.t_adjusts = tmp(:);
    tmp = load('timing_files_201/t_delays.txt');
    timing_parms.t_delay = tmp(:);
    %}
    timing_parms.t_tag                  = 0.05 *ones(Nframes,1) ;
    timing_parms.ArtSup_delay           = 0.100 *ones(Nframes,1) ;
    timing_parms.t_aq                   = 0.0329*18*ones(Nframes,1) ;
    timing_parms.order                  = 1;
    
    delays = timing_parms.t_delay;
    ordr = timing_parms.order;
    
    parms.mtis0 =     1 ;
    parms.Disp =      40;
    parms.r1blood = 1/1.7;
    
    parms.f =       0.01;
    parms.cbva =    0.01;
    parms.kfor =    0;
    parms.bat2 =    0;
    parms.bat =     0.05;
    parms.r1tis =   1/0.9;%1.4;
    parms.flip =    deg2rad(90);
    parms.r2tis=    1/0.090;
    
    % generate Data : single run
    test_signal = abs(single(  ...
        gen_signals_vs_211024(parms,...
        delays, ...
        timing_parms,...
        doFigs,doSub) ...
        + (0*(randn(1,Nframes))))); % LH Model one BAT
    
    parms.f = 0.005;
    test_signal2 = abs(single(  ...
        gen_signals_vs_211024(parms,...
        delays, ...
        timing_parms,...
        doFigs,doSub) ...
        + (0*(randn(1,Nframes))))); % LH Model one BAT
    
    %{
    plot(test_signal)
    hold on
    plot(test_signal2)
    title(sprintf('RMS change from CBF increase (50%%) : %0.2e',...
        norm((test_signal-test_signal2)./test_signal)));
    %    norm((test_signal-test_signal2))./norm(test_signal)));
    
    hold off
    drawnow
    %}
    rms_change = norm((test_signal-test_signal2))/norm(test_signal);
    if rms_change > best
        rp_best= rp;
        rp2_best = rp2;
        best = rms_change;
    end
end

% now test the signal wit the best RMS
timing_parms.labelcontrol       = ones(Nframes,1);
timing_parms.labelcontrol(rp_best)   = 0;
timing_parms.doArtSup           = ones(Nframes,1);
timing_parms.doArtSup(rp2_best)      = 0;

% generate Data : single run
parms.f = 0.01;

test_signal = abs(single(  ...
    gen_signals_vs_211024(parms,...
    delays, ...
    timing_parms,...
    doFigs,doSub) ...
    + (0*(randn(1,Nframes))))); % LH Model one BAT

parms.f = 0.005;
test_signal2 = abs(single(  ...
    gen_signals_vs_211024(parms,...
    delays, ...
    timing_parms,...
    doFigs,doSub) ...
    + (0*(randn(1,Nframes))))); % LH Model one BAT

subplot(311)
plot(timing_parms.t_delay); title('t delay')
subplot(312)
plot(timing_parms.t_adjusts); title('t adjust')

subplot(313)
plot(test_signal)
hold on
plot(test_signal2)
plot(100*abs(test_signal - test_signal2));

axis([1 length(test_signal) 0 1])
title(sprintf('NRMS change from 50%% CBF increase : %0.2e',...
    best));
hold off


seqDuration = sum(timing_parms.t_delay)+ ...
    sum(timing_parms.t_tag) + ...
    sum(timing_parms.t_adjusts) + ...
    sum(timing_parms.ArtSup_delay)+...
    sum(timing_parms.t_aq)

text(10, double(min(test_signal2)), sprintf('Duration= %0.2f seconds', seqDuration));
drawnow

%%

f       = rand(1,nData) * 0.025;
cbv     = rand(1,nData) * 0.05;
kfor    = rand(1,nData) * 0;
bat_tis = rand(1,nData) * 0;
bat_art = rand(1,nData) * 0.5 + 0.01;
r1tis   = 0.3 + 2*rand(1,nData) ;
flip    = (pi/2)*( 0.9 + 0.2*rand(1,nData));
r2tis   = 9 + (16*rand(1,nData));  % T2 values from 40 ms to 111ms 


theta = [
    f;
    cbv;
    kfor;
    bat_tis;
    bat_art;
    r1tis;
    flip;
    r2tis;
    
    ]; % LH model

theta = single(theta)';


delays = timing_parms.t_delay;


%%
fprintf('\nBegin generating %d time courses', nData);
%parpool('local', 'SpmdEnabled', false)
%pp =  parpool(24)
doFigs = 0;

parfor idx = 1:nData
    
    %%%% flip this for
    parms = arr2parms(theta(idx,:));
    
    Y(idx,:) = abs(single(  ...
        gen_signals_vs_211024(parms,...
        delays,...
        timing_parms,...
        doFigs,doSub) ...
        + (sigma*(randn(1,Nframes))))); % LH Model one BAT
    
    % plot(Y(idx,:))
    % drawnow

    if mod(idx,1e5)==0
        fprintf('\ngenerating time course number ... %e of %e', idx , nData);
    end
end

fprintf('\n ....Finished generating %d time courses', nData);



%% split into training and testing data

lim = round(split*nData);
trainData = Y(1:lim,:);
trainTgts = theta(1:lim,:);

testData = Y(lim+1:end,:);
testTgts = theta(lim+1:end,:);
%%

!mkdir chirpfun_vss/
save -v7.3 chirpfun_vss/trainData.mat trainData
save -v7.3 chirpfun_vss/trainTgts.mat trainTgts
save -v7.3 chirpfun_vss/testData.mat testData
save -v7.3 chirpfun_vss/testTgts.mat testTgts

return

function parms = arr2parms(theta)

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       theta(1);
parms.cbva =    theta(2);
parms.kfor =    theta(3);
parms.bat2 =    theta(4);
parms.bat =     theta(5);
parms.r1tis =   theta(6);
parms.flip =    theta(7);

% R2 of tissue.  WIll determine the efficiency/effect on tissue and arteries
parms.r2tis = theta(8);


return
