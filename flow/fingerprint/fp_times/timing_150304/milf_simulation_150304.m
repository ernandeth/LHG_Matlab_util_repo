Nframes = 80;
nslices = 1;
psdseqtime = 0.030;
% read in sequence timings
%{
label_dur = [linspace(0.2, 2, 20) linspace(2,2,10) linspace(2,0.2, 20) ];
label_dur(2:2:end) = label_dur(1:2:end);
save t_tags.txt label_dur -ascii

PID =[linspace(0.05, 0.05,20) linspace(0.05, 2, 10) linspace(0.05, 1,20) ];
PID(2:2:end) = PID(1:2:end);

save t_delays.txt PID -ascii

t_adjust = [linspace(1, 0.2, 20)  linspace(2, 0.05, 10)    linspace(0.05, 0.05, 20) ];
t_adjust(2:2:end) = t_adjust(1:2:end);
%%
%}

label_dur = [linspace(0.2, 2, 20) linspace(0.2, 2, 20)];
%label_dur(2:2:end) = label_dur(1:2:end);
label_dur = repmat(label_dur,1,2);
save t_tags.txt label_dur -ascii

PID =[linspace(2.0, 0.05,20)   linspace(0.05, 2,20) ];
%PID(2:2:end) = PID(1:2:end);
PID = repmat(PID,1,2);

save t_delays.txt PID -ascii

t_adjust = [linspace(0.1, 1, 20) linspace(1.0, 0.05, 20) ];
%t_adjust(2:2:end) = t_adjust(1:2:end);
t_adjust = repmat(t_adjust,1,2);

t_adjust = 4 - PID - label_dur;

save t_adjusts.txt t_adjust -ascii

Nframes = length(t_adjust);


% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to
PID = PID + 0.010;
t_adjust = t_adjust + 0.008;
label_dur = label_dur(1:Nframes);
t_adjust =  t_adjust(1:Nframes);
PID = PID(1:Nframes);

TR = t_adjust +  label_dur + PID + nslices*psdseqtime ;
%
duration = sum(TR)+2;

timing_parms.PID = PID;
timing_parms.label_dur = label_dur;
timing_parms.t_adjust = t_adjust;
timing_parms.order = 1; % tag-control
timing_parms.order = 2; % control-tag

%%
% set up default values for fit
phys_parms.f =         0.008;
phys_parms.mtis0 =     1 ;
phys_parms.cbva =       0.01
phys_parms.transit=    1.5 ;
phys_parms.kfor =     0.02 ;
phys_parms.r1tis =     1/1.4  ;
phys_parms.beta =      70*pi/180 ; % flip angle in radians
phys_parms.L = 1;
phys_parms.Disp =      25;
phys_parms.Ptime =     0.5;

obs = gen_signals_140618(phys_parms, timing_parms, 1, 0);
[dict, parms] = gen_dictionary_140618(timing_parms);
xc = (corrcoef(dict'))

xc = (xc(:));
G = mean(xc)