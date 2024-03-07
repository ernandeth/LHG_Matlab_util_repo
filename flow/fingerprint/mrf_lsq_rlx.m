function [r1 r2 b1err] = mrf_lsq_rlx(timeseries, aq_parms, Fparms, showFit)
% function [r1 r2 b1err] = mrf_lsq_rlx(timeseries, aq_parms, Fparms, showFit)
%
% Estimate relaxation parameters from VSASL MRF data
% Do a least squares fit of the data to a model of VSASL
% MRF data for a specific set of aquisition parameters (schedule)
%
% It uses the gen_signals_vs_230321 for the model
% New: It uses the gen_signals_vs_230718 for the model
%
% INputs: 
%   timeseries: the data to fit
%   aq_parms:   structure with readout parameters
%   Fparms :    flow related parameters: cbf cbv, bat
%   showFit:    this parameter lets you watch the fit evolve.  It's for
%               troubleshooting only
%
% It estimates only R1 and R2 and M0,  assuming the other parameters
% in the model
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
%  [r1 r2 b1err] = mrf_lsq(timeseries, aqparms, 0)
%
% For a demo, you can also call it without any parameters
%

if nargin==0
    Fparms = [0.01 0.02 0.2];
end

% Set some Tissue parameters
parms.f=          Fparms(1);
parms.cbva =      Fparms(2) ;
parms.bat =       Fparms(3) ;
parms.Mtis0 =     1 ;
parms.r1tis =     0.7  ;
parms.r2tis =     11 ;
parms.b1err =   0.0;

if nargin==0
    % Set pulse sequence parameters
    Nframes = 20;
    aq_parms.label_type =    'BIR8'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    aq_parms.RO_type =      'GRE'; %'FSE';   % 'GRE'
    aq_parms.flip =         deg2rad(30);
    aq_parms.t_tags =        0;% 0.1*ones(Nframes,1);
    aq_parms.del1 =          0.2*ones(Nframes, 1);
    aq_parms.del2 =          0.2+0.5*sin(linspace(0.1,4*pi,Nframes))';
    aq_parms.del3 =          0.1*ones(Nframes,1);  % delay between AS pulse and acqusition
    aq_parms.labelcontrol =  zeros(Nframes,1);
    aq_parms.labelcontrol(2:2:end)= 1;
    aq_parms.order =         1;
    aq_parms.doArtSup =      ones(Nframes,1);
    aq_parms.doArtSup(1)=    0;
    aq_parms.t_aq = 0.750 ;  % duration of the whole readout

    aq_parms.Ma = [];

    % generate test data
    doSub = 0;
    dofigs = 0;
    %data = gen_signals_vs_230718(parms, aq_parms, dofigs,doSub);
    [data Mart]= gen_signals_vs_230918(parms, aq_parms, dofigs,doSub);
    aq_parms.Ma = Mart;
    
    timeseries = data + 0.01*randn(size(data));
    showFit = 1;

end


% checking the model fit for this schedule using synth data
% (if you uncomment this, things will fit really well)
%{
%data = gen_signals_vs_230321(parms, aq_parms, 0,0);
data = gen_signals_vs_230718(parms, aq_parms, 0,0);
timeseries = data + 0.01*randn(size(data));
showFit=1
%}

% How to scale the data??
timeseries = timeseries/norm(timeseries);

% set up the estimation problem
% typical R1 , R2,  B1error values
guess0= [1/1.2  1/0.08  0.05];
lb =    [1/4    1/0.5   -0.25];
ub =    [1/0.4  1/0.02  0.25]; 

opts = optimoptions(@lsqnonlin);
%opts.MaxIterations = 1000;
opts.SubproblemAlgorithm="cg";
% opts.FiniteDifferenceStepSize = 1e-8;
% %opts.MaxFunctionEvaluations = 1e3;
% opts.OptimalityTolerance = 1e-8;
% opts.FunctionTolerance = 1e-8;
% opts.StepTolerance = 1e-8;
opts.FiniteDifferenceType ='central';
opts.Display='off';

if showFit
    opts.Display = "iter-detailed";
end
%
[est res] = lsqnonlin(@(guess)myfun...
    (guess, aq_parms, parms, showFit, timeseries), ...
    guess0, ...
    lb,ub, opts);
%{
[est res] = fsolve(@(guess)myfun(...
    guess, aq_parms, parms, showFit, timeseries), ...
    guess0);
%}
r1 = est(1);
r2 = est(2);
b1err = est(3);

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

r1tis = guess(1);
r2tis = guess(2);
b1err = guess(3);

parms.r1tis = r1tis;
parms.r2tis = r2tis;
parms.b1err = b1err;


%synthdata = gen_signals_vs_230718(parms, aqparms, 0, 0);
synthdata = gen_signals_vs_230918(parms, aqparms, 0, 0);

% How to scale the data??
synthdata = synthdata/norm(synthdata);

% res =(data-synthdata).^2  ;  %  <--- good for synthetic data
res = data-synthdata  ;

if showFit
    plot(data,'*')
    hold on
    plot(synthdata)
    hold off
    legend('Data', 'Guess')
    title(sprintf('R1 %0.2f, R2 %0.2f B1err %d', parms.r1tis, parms.r2tis, parms.b1err))
    
    drawnow
end
return