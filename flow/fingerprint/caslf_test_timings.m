psdseqtime = 0.03;
nslices = 1;

label_dur = load('t_tags.txt');
PID = load('t_delays.txt');
t_adjust = load('t_adjusts.txt');

% corrections - pulse sequence in the scope doesn't quite do what it's
% supposed to
PID = PID + 0.010;
t_adjust = t_adjust + 0.008;

timing_parms.PID = PID;
timing_parms.label_dur = label_dur;
timing_parms.t_adjust = t_adjust;

[dy p]= gen_dictionary_test(timing_parms);

xc = (corrcoef(dy'))

plot(dy')
