% Generates ASL Training Data for a machine learning based estimation
% technique. The parameters are uniformly random over a range,
%
% evolved from AL code for VSI training
% LHG: jan 2022
% modiifcations for new VSASL model.  Key difference:  reduced number of
% parameters.  (No kfor, bat2 ... alphas)
%
% All the labeling efficiencies are calculated from T2 of the
% blood and tissue, and the type of pulse.
% - Note that each  type VS pulse has an effective TE

useFileSchedule = 0;

% number of training + testing time-series and the length of time series;
nData = 6.0*10^6;
nData = 10^6;

% nData = 50;  % small data set to make sure things work right


% make up our own schedule on the fly
%     load optimal_sequence.mat
%     timing_parms = best_timing_parms;
load timing_parms.mat
Nframes                     = length(timing_parms.t_delay);


% synthetic training data matrix
Y = single(zeros(nData,Nframes));

% train-test split ratio
split = 0.9;

%
delays = timing_parms.t_delay;

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

parms.mtis0 =   1 ;
parms.Disp =    40;
parms.r1blood = 1/1.7;

parms.f =       0.005;
parms.cbva =    0.005;
parms.bat =     0.01;
parms.r1tis =   1/1.2;
parms.r2tis =   1/0.050;
parms.flip =    deg2rad(30);

test_signal = abs(single(  ...
    gen_signals_vs_211024(parms,...
    delays, ...
    timing_parms,...
    doFigs,doSub)));


seqDuration = sum(timing_parms.t_delay)+ ...
    sum(timing_parms.t_tag) + ...
    sum(timing_parms.t_adjusts) + ...
    sum(timing_parms.ArtSup_delay)+...
    sum(timing_parms.t_aq)

plot(test_signal)
title(sprintf('Sequence Duration: %f seconds', seqDuration));
drawnow
%%
% Now generate distribution of tissue parameters
%{
order:
    parms.f =       p_array(1);
    parms.cbva =    p_array(2);
    parms.bat =     p_array(3);
    parms.r1tis =   p_array(4);
    parms.r2tis =   p_array(5);
    parms.flip =    p_array(6);

%}
f       = rand(1,nData) * 0.02;
cbv     = rand(1,nData) * 0.05;
bat_tis = 0.01 + rand(1,nData) * 1;
r1tis   = 0.3 + 2*rand(1,nData) ;
r2tis   =  10 + 40*rand(1,nData) ;
flip    =  pi/20 + (pi/2)*rand(1,nData);

theta = [
    f;
    cbv;
    bat_tis;
    r1tis;
    r2tis;
    flip];

theta = single(theta)';


delays = timing_parms.t_delay;



%%
fprintf('\nBegin generating %d time courses', nData);
%parpool('local', 'SpmdEnabled', false)
%pp =  parpool(24)
doFigs = 0;

parfor idx = 1:nData
    
    parms = make_parm_struct(theta(idx,:));
    
    Y(idx,:) = abs(single(  ...
        gen_signals_vs_211024(...
        parms,...
        delays,...
        timing_parms,...
        doFigs,doSub)));
    
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
%{
!mkdir AL_schedule_model_220130/
save -v7.3 AL_schedule_model_220130/trainData_2.mat trainData
save -v7.3 AL_schedule_model_220130/trainTgts_2.mat trainTgts
save -v7.3 AL_schedule_model_220130/testData_2.mat testData
save -v7.3 AL_schedule_model_220130/testTgts_2.mat testTgts
toc
%}
!mkdir chirpfun_vss2
save -v7.3 chirpfun_vss2/trainData.mat trainData
save -v7.3 chirpfun_vss2/trainTgts.mat trainTgts
save -v7.3 chirpfun_vss2/testData.mat testData
save -v7.3 chirpfun_vss2/testTgts.mat testTgts
save chirpfun_vss2/timing_parms.mat timing_parms
toc
%%
function parms = make_parm_struct(p_array)

parms.f =           p_array(1);
parms.cbva =        p_array(2);
parms.bat =         p_array(3);
parms.r1tis =       p_array(4);
parms.r2tis =       p_array(5);
parms.flip =        p_array(6);
parms.mtis0 =       1 ;
parms.Disp =        40;
parms.r1blood =     1/1.7;


end
