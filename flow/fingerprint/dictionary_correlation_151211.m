function G =dictionary_correlation_151211 (x)
% function G = dictionary_correlation_151211(x)
%
% this function  takes in a vector of timing parameters and uses it to
% create a dictionary of ASL signals.
% it returns the correlation among the dictionary entries
%
%
% it returns the norm of the correlation matrix of the dictionary
%
% maxTR = 100/length(x); % 100 second experiment.

global isLabel_global

maxTR = 100/length(x);  % keep this at 100 seconds, so 0.5 seconds per TR

tvec = x;

N = length(tvec);

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

%
% Combinations of  Physiological parameters for parameters:
dict_phys_parms.r1tis = 1/1.5; % 1./linspace(0.5, 3, 3);
dict_phys_parms.flip =   deg2rad(60); % deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

% For debugging:
if 0
    phys.r1tis = 1;
    phys.flip =  deg2rad(60);
    phys.bat =  1;
    phys.f =  0.01;
    phys.cbva = 0.01;
    phys.kfor = 0.02;
    phys.Disp = 20;
    phys.mtis0 = 1;
    entry = gen_signals_150521( phys , timing_parms, 1, 0);
    %plot(timing_parms.t_tag); drawnow
end

[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

sigma2 = sum(var(dict'));

xc = abs(corrcoef(dict'));

% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

% penalize entries greater than 0.9 correlation
xc(xc>0.9) = 1.5*xc(xc>0.9);

xc = (xc(:));
% G = mean(xc);
% I'm adding a penalty for having too little variance.
G = mean(xc) + 0.5/sigma2;
return
    
    
