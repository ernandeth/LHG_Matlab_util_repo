% 
% ObjectiveFunction = @dictionary_correlation_151211
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjectiveFunction = @dictionary_correlation_151211;
nvars = 200;    % Number of acquisitions - 
nvars = 100;
% maxTR = 2*0.56;
maxTR = 100/nvars;  % keep this at 100 seconds, so 4 seconds per TR
% maxTR = 20/nvars;  % keep this at 20 seconds, so 0.5 seconds per TR
global isLabel_global

initguess = 0.89 * rand(500,nvars) + 0.01;

% constrain the label durations:
% lbounds(1:nvars) = 0;
% ubounds(1:nvars) = 2.5;   
% ubounds(1:nvars) = 0.8;   

lbounds = 0.01 * ones(1,nvars);
ubounds = 0.9 * ones(1,nvars);

opts = gaoptimset;

opts.Generations = 100;
opts.PopulationSize = 10*nvars;
opts.PlotFcns = @gaplotbestf;
opts.Display = 'iter';
opts.StallGenLimit = 100;
% opts.PopulationType = 'bitstring';
options.TolCon=1e-16;
options.EliteCount = 5; %nvars/20;
opts.InitialPopulation = initguess;
options.CrossoverFraction=.9; %default.8
options.MutationFcn={@mutationgaussian 1 1};
%options.MutationFcn={@mutationadaptfeasible};
options.UseParallel= true;

% initguess = ones(1,nvars);
% initguess(2:2:end) = 0;



% force the t_tag to some initial randomized order 
Nframes = nvars;
isLabel_global = round(rand(1, Nframes));


% Call the optimizer:
[x,fval] = ga(ObjectiveFunction, nvars, [],[], [], [], lbounds, ubounds,[],[],opts);

tvec = x;

N = length(tvec);

% code from the objectve function here:
%{
timing_parms.t_tag =  0.4 * ones(size(tvec));
timing_parms.t_adjust = 0.01* ones(size(tvec)) ;
timing_parms.t_adjust =0.05 * ones(1,N) ;
timing_parms.isLabel = x;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms. t_aq =ones(1,N) * 0.035;
timing_parms.t_delay = maxTR -timing_parms.t_tag - timing_parms.t_adjust - timing_parms.t_aq;  
%}

% optimizing sequence of t_tags
timing_parms.t_tag = tvec;
% using preshuffled order of controls and labels.
timing_parms.isLabel = isLabel_global;



timing_parms.t_delay = ones(size(tvec)) * 0.01;
timing_parms.t_aq =ones(size(tvec)) * 0.035;
timing_parms.t_adjust = maxTR - timing_parms.t_tag - timing_parms.t_delay - timing_parms.t_aq;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;

timing_parms.order = 1;
timing_parms.Nlabel_group = 1;


TR = timing_parms(:).t_adjust +...
    timing_parms(:).t_tag +...
    timing_parms(:).t_delay +...
    timing_parms(:). t_aq ;
ExperimentLength = sum(TR)

% Combinations of  Physiological parameters for parameters:
dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

figure

%%

% Now show me what this dictionary does
[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

plot(dict')
drawnow

xc = (corrcoef(dict'));
xc = abs(xc(:));
G = mean(xc)
G = norm(xc);

save optimalTiming_parms dict parms timing_parms G x

figure
subplot(211)
imagesc(dict)

xc = (corrcoef(dict'));
subplot(212)
imagesc((xc))

xc = (xc(:));
G = mean(xc)

%%
% Luis Gomez's options:
%{
options=gaoptimset;

options.CrossoverFraction=.6;%default.8
options.PopInitRange(1,:)=-999*ones([1 size]);%options.PopInitRange=[0;1]
options.PopInitRange(2,:)=999*ones([1 size]);
options.PopulationSize=5000;%options.PopulationSize=3000
options.PopulationType='doubleVector';%'bitString';%'doubleVector';

xsymm=5000*rand([options.PopulationSize size]);

options.EliteCount=5;%options.EliteCount=2
options.MigrationDirection='forward';%options.MigrationDirection='forward'
options.MigrationInterval=10;%options.MigrationInterval=20
options.MigrationFraction=0.1;%options.MigrationFraction=0.2
options.Generations=100;%options.Generations=200
options.TimeLimit=inf;%options.TimeLimit=inf
options.FitnessLimit=-inf;%options.FitnessLimit=-inf
options.StallGenLimit=200000;%options.StallGenLimit=1000
options.StallTimeLimit=20000;%options.StallTimeLimit=20
options.TolFun=1e-6;%options.TolFun=1e-6
options.TolCon=1e-6;%options.TolCon=1e-6
options.InitialPopulation=xsymm;%xsymm;%options.InitialPopulation=[]
options.InitialScores=[];%options.InitialScores=[]
options.InitialPenalty=10;%options.InitialPenalty=10
options.PenaltyFactor=100;%options.PenaltyFactor=100
options.PlotInterval=10;%options.PlotInterval=1
options.CreationFcn=@gacreationuniform;%options.CreationFcn=@gacreationuniform
options.FitnessScalingFcn=@fitscalingrank;%options.FitnessScalingFcn=@fitscalingrank

options.SelectionFcn=@selectionstochunif
options.CrossoverFcn=@crossoverscattered;%arithmetic;%options.CrossoverFcn=@crossoverscattered
options.MutationFcn={@mutationgaussian 1 1};%options.MutationFcn={{@mutationgaussian 1  1}}
options.HybridFcn=[];%options.HybridFcn=[]
options.Display='iter';%options.Display='final'
options.PlotFcns=[];%options.PlotFcns=[]
options.OutputFcns={[@gaoutputgen1]};
options.Vectorized='off';%options.Vectorized='off'
%}
