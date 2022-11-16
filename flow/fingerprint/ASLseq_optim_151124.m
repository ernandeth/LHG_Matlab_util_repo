% 
% ObjectiveFunction = @dictionary_correlation_150918
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ObjectiveFunction = @dictionary_correlation_150918;
nvars = 50;    % Number of acquisitions - 2 variables at each iteration
maxTR = 3;

% constrain the label durations:
lbounds(1:nvars) = 0;
ubounds(1:nvars) = 1;   

% constrain the t_delay s
lbounds(nvars+1 : 2*nvars) = 0.05;
ubounds(nvars+1 : 2*nvars) =1;

%
initguess =  rand(1,nvars*2);


% Call the optimizer:
[x,fval] = lsqnonlin(ObjectiveFunction, initguess); %, lbounds, ubounds);

tvec = x;

N = length(tvec)/2;
timing_parms.t_tag =  tvec(1 : N);
timing_parms.t_delay = tvec(N+1 : 2*N) ;

timing_parms.t_adjust =0.05 * ones(1,N) ;
timing_parms.isLabel = ones(N,1)  ;
timing_parms.isLabel(2:2:end) = 0  ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms. t_aq =ones(1,N) * 0.035;
%timing_parms.t_adjust = maxTR -timing_parms.t_tag - timing_parms.t_delay - timing_parms.t_aq;  
TR = timing_parms(:).t_adjust +...
    timing_parms(:).t_tag +...
    timing_parms(:).t_delay +...
    timing_parms(:). t_aq ;
ExperimentLength = sum(TR)

% Combinations of  Physiological parameters for parameters:
dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);
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
