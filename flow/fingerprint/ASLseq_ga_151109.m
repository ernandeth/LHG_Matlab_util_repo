% 
% ObjectiveFunction = @dictionary_correlation_151109
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjectiveFunction = @dictionary_correlation_151109;
nvars = 40;    % Number of acquisitions - 2 variables at each iteration
maxTR = 2*0.56;
maxTR = 100/nvars;  % keep this at 100 seconds, so 4 seconds per TR
maxTR = 20/nvars;  % keep this at 20 seconds, so 0.5 seconds per TR

% constrain the label durations:
lbounds(1:nvars) = 0;
ubounds(1:nvars) = 2.5;   
ubounds(1:nvars) = 0.4;   

% constrain the t_delay s
%lbounds(nvars+1 : 2*nvars) = 0.05;
%ubounds(nvars+1 : 2*nvars) = 0.4;

opts = gaoptimset;

opts.Generations = 10*nvars;
opts.PopulationSize = 25*nvars;
opts.PlotFcns = @gaplotbestf;
opts.Display = 'iter';
opts.StallGenLimit = 5;

%
% initguess =  1.8*rand([opts.PopulationSize  nvars]);
initguess =  0.4*rand([opts.PopulationSize  nvars]);

options.EliteCount = 6;
opts.InitialPopulation =initguess;


% Call the optimizer:
[x,fval] = ga(ObjectiveFunction, nvars, [],[], [], [], lbounds, ubounds,[],[],opts);

tvec = x;

N = length(tvec);

timing_parms.t_tag =  tvec(1 : N);
timing_parms.t_adjust = 0.01* ones(size(tvec)) ;

timing_parms.t_adjust =0.01 * ones(1,N) ;

timing_parms.isLabel = ones(1, N)  ;
timing_parms.isLabel(2:2:end) = 0  ;

%timing_parms.isLabel = ones(1, N)  ;
%timing_parms.isLabel(5:8) = 0  ;
%timing_parms.isLabel(13:16) = 0  ;

timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms. t_aq =ones(1,N) * 0.035;
timing_parms.t_delay = maxTR -timing_parms.t_tag - timing_parms.t_adjust - timing_parms.t_aq;  



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
