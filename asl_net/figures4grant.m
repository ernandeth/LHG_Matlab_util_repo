addpath  /home/hernan/matlab/flow/fingerprint

Nframes = 50


t_tag = ones(Nframes,1 ) * 1.8;
t_delay =  ones(size(t_tag)) * 1.8;
isLabel = ones(size(t_tag));
isLabel(1:2:end) = 0;
isLabel(1:8) = 0;
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = ones(size(t_tag)) * 4  ;
t_adjust = t_adjust - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

parms0 = struct( ...
    'mtis0', 1,...
    'f', 0.01 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0, ...
    'r1tis', 0.8,...
    'flip', 0, ... %pi/4, ...
    'Disp', 20);


% show the signal for one case
doSub = 0;
dofigs = 1;

entry = gen_signals_180320(parms0, timing_parms, dofigs, doSub);

% in the breakpont before dofigs
plot(timevec, Ma)
axis([50 60 0 1])
ylabel('M_z')
xlabel('Time (sec)')
title('Arterial Magnetization')
fatline1

print -dpng Ma

plot(timevec, M)
axis([50 60 0.98 1])
ylabel('M_z')
xlabel('Time (sec)')
title('Tissue Magnetization')
fatlines
print -dpng Mt


Ma = labelfun;
Ma(Ma<0.1) = -1;
plot(timevec, -Ma*0.9)
axis([50 60 -1 1])
ylabel('M_z')
xlabel('Time (sec)')
title('Carotid Artery Magnetization')
fatlines
print -dpng Ma_carotid

%%
