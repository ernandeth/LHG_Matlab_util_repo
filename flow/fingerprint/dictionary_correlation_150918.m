function G =dictionary_correlation_151108 (x)
% function G = dictionary_correlation_151108(x)
%
% this function  takes in a vector of timing parameters and uses it to
% create a dictionary of ASL signals.
% it returns the correlation among the dictionary entries
%
%
% it returns the norm of the correlation matrix of the dictionary
%
maxTR = 3;
tvec = x;

N = length(tvec);
timing_parms.t_tag =  tvec(1 : N);
timing_parms.t_delay = 0.01* ones(size(tvec)) ;

timing_parms.t_adjust =0.05 * ones(1,N) ;
timing_parms.isLabel = ones(N,1)  ;
timing_parms.isLabel(2:2:end) = 0  ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms. t_aq =ones(1,N) * 0.035;
timing_parms.t_adjust = maxTR -timing_parms.t_tag - timing_parms.t_delay - timing_parms.t_aq;  

%
% Combinations of  Physiological parameters for parameters:
dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);


xc = abs(corrcoef(dict'));
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

% penalize entries greater than 0.9 correlation
xc(xc>0.9) = 1.5*xc(xc>0.9);

xc = (xc(:));
G = mean(xc);
%G = norm(xc);
return
    
    
