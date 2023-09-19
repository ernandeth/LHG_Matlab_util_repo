function [est, R] = mrf_lsq_flow(timeseries, aq_parms, Rparms, showFit)
% function [ [cbf cbv bat], residual] = mrf_lsq_flow(timeseries, aq_parms, Rparms, showFit)
%
% Estimates FLow parameters from VSASL MRF data
% Do a least squares fit of the data to a model of VSASL
% MRF data for a specific set of aquisition parameters (schedule)
%
% It uses the gen_signals_vs_230718 for the model
% New: It uses the gen_signals_vs_230321 for the model
%
% INputs: 
%   timeseries: the data to fit
%   aq_parms:   structure with readout parameters
%   Rparms:     structure with precomputed parameters - usually the pre-estimated R1 and
%               R2 from a previous estimation step:
%                   [r1 r2 ]
%   showFit:    this parameter lets you watch the fit evolve.  It's for
%               troubleshooting only
%
% It estimates only CBF and CBV
%
% example:
%  cd MRF_data_dir
%  aqparms = read_timing_files('./')
%
%  raw = readnii('timeseries_mag.nii');
%  msk = makevarmask(raw , 50);
%  raw = raw.*msk;
%  timeseries = squeeze(sum(raw, [1,2,3]));
%
%  [r1 r2] = mrf_lsq(timeseries, aqparms, 0)
%
%  [cbf cbv], res = mrf_lsq_flow(timeseries, aqparms, parms 0)
%
% For a demo, you can also call it without any parameters
%

% Set some Tissue parameters
parms.f=          0.015;
parms.cbva =      0.02 ;
parms.bat =       0.15 ;
parms.Mtis0 =     1 ;
parms.flip =      deg2rad(90) ; % flip angle in radians
parms.r1tis =     0.7  ;
parms.r2tis =     11 ;
parms.b1err =       0;


if nargin==0
    % Set pulse sequence parameters
    Nframes = 20;
    aq_parms.label_type =    'BIR8inv'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    aq_parms.RO_type =       'FSE'; %'FSE';   % 'GRE'
    aq_parms.t_tags =        0;% 0.1*ones(Nframes,1);
    aq_parms.del1 =          0.5*ones(Nframes, 1);
    aq_parms.del2 =          0.2 + 1.2*abs(sin(linspace(0.1,pi*2,Nframes)))';
    aq_parms.del3 =          0.1*ones(Nframes,1);  % delay between AS pulse and acqusition
    aq_parms.labelcontrol =  zeros(Nframes,1);
    aq_parms.labelcontrol(2:2:end)= 1;
    aq_parms.order =         1;
    aq_parms.doArtSup =      ones(Nframes,1);
    aq_parms.doArtSup(1)=    0;
    aq_parms.t_aq =          0.750 ;  % duration of the whole readout

aq_parms = read_timing_files('./')
    
% Relaxation parms from user input
    Rparms(1)= parms.r1tis;
    Rparms(2) = parms.r2tis;
    Rparms(3) = parms.b1err;

    if(aq_parms.RO_type=='GRE')
        parms.flip = deg2rad(20);
    end

    % generate test data
    doSub = 0;
    dofigs = 1;
    %data = gen_signals_vs_230718(parms, aq_parms, dofigs,doSub);
    data = gen_signals_vs_230918(parms, aq_parms, dofigs,doSub);
    
    timeseries = data + 0.002*randn(size(data)) + 0.01*sin(linspace(-2,4,length(data)));
    showFit = 1;

end

% use the predefined relaxation times here
parms.r1tis = Rparms(1);
parms.r2tis = Rparms(2) ;
parms.b1err = Rparms(3);

% use this for GRE
%parms.flip = deg2rad(10);
%aq_parms.RO_type =        'GRE'

% checking the model fit for this schedule using synth data
% (if you uncomment this, things will fit really well)
%{
data = gen_signals_vs_230321(parms, aq_parms, 0,0);
timeseries = data + 0.01*randn(size(data));
%}

% How to scale the data??
timeseries = timeseries/norm(timeseries);

% set up the estimation problem
guess0= [0.01  0.01  0.05];%  typical cbf, cbva, bat 
lb =    [0.001  0.0   0.01];
ub =    [0.05   0.1   0.5]; 

opts = optimoptions(@lsqnonlin);

%opts.MaxIterations = 1000;
opts.SubproblemAlgorithm="cg";
% opts.FiniteDifferenceStepSize = 1e-8;
% %opts.MaxFunctionEvaluations = 1e3;
% opts.OptimalityTolerance = 1e-8;
% opts.FunctionTolerance = 1e-8;
% opts.StepTolerance = 1e-8;
opts.FiniteDifferenceType ='central';

if showFit
    opts.Display = "iter-detailed";
end
%{
[est R] = lsqnonlin(@(guess)myfun...
    (guess, aq_parms, parms, showFit, timeseries), ...
    guess0, ...
    lb,ub, opts);
%}
[est R] = fmincon(@(guess)myfun...
    (guess, aq_parms, parms, showFit, timeseries), ...
    guess0, [],[],[],[],...
    lb,ub);
cbf = est(1);
cbv = est(2);
bat = est(3);


return

%%
function res = myfun(guess, aqparms, parms, showFit, data )
% function res = myfun(guess, aqparms, parms, data )
%
% returns the difference between the input data
% and the data generated by a model , whose inputs are specified
% by guess, aqparms and parms
%
% this function is meant to be called by lsqnonlin
%

%% Testing: 8.10.23
doSub = 0;
if doSub
    data = data(1:2:end)- data(2:2:end);
end
%%
cbf = guess(1);
cbv = guess(2);
bat = guess(3);

parms.f = cbf;
parms.cbva = cbv;
parms.bat = bat;

%synthdata = gen_signals_vs_230718(parms, aqparms, 0, doSub);
synthdata = gen_signals_vs_230918(parms, aqparms, 0, doSub);

% How to scale the data??
synthdata = synthdata/norm(synthdata);

%res =(data-synthdata).^2  ;  %  <--- good for synthetic data

% in fmincon you can specify the exact cost function
res = norm(data-synthdata)^2  +  1e-5*norm(data-synthdata,1);

% res = diff(res);  %<---- very brute force highpass filter
% res = data .* synthdata';  % looking for a dictionary match by pattern
% recognition

if showFit
    plot(data,'*')
    hold on
    plot(synthdata)
    hold off
    legend('Data', 'Guess')
    title(sprintf('R1 %0.2f, R2 %0.2f M0 %d', parms.r1tis, parms.r2tis, parms.Mtis0)) 
    drawnow
end

return