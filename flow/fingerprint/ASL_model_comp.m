close all; clear variables;

% Sequence Parameters for Testing
parms.mtis0=1;
parms.f=50/6000;
parms.cbva=0.0001;
parms.transit=1.2;
parms.kfor=0;
parms.r1tis=1/1.4;
parms.beta=30;
parms.Disp=10;

PIDval = 1.6;
ldval = 3;
alpha = 0.9;
lpc = 0.9;
T1blood = 1.67;
obsrep = 10;
timing_parms.PID = ones(1,obsrep).*PIDval;
timing_parms.label_dur = ones(1,obsrep).*ldval;
timing_parms.imTR = 0.03;
timing_parms.TR = timing_parms.PID+timing_parms.label_dur+timing_parms.imTR;
timing_parms.imTR = 1;
timing_parms.TE = 2.5/1000;

obs = gen_signals_140618_kwv2(parms, timing_parms, 1, 0);
diff_m=(obs(end-1)-obs(end))./sin(parms.beta)

t1=ldval+PIDval;
R1p = parms.r1tis + parms.f/lpc;
T1p = 1/R1p;
diff_SI = 2*parms.mtis0*parms.f*T1p*alpha*exp(-parms.transit/T1blood)*exp(-(t1-ldval-parms.transit)*R1p)*...
    (1-exp(-(ldval)*R1p))

diff_SI/diff_m
