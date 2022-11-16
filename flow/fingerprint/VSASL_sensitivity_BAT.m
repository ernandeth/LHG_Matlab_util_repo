%% creating ASL timing schedule for BAT measurement
t_aq        = 0.51;

Nframes = 50;


isLabel     = ones(Nframes,1);
isLabel(2:2:end) = 0;

timing_parms.doArtSup       = 1 *ones(Nframes,1);
timing_parms.ArtSup_delay   = 0.05 * ones(Nframes,1);
timing_parms.Nlabel_group   = 1;
timing_parms.t_tag          = 0.3 *ones(Nframes,1 );

t_delay  = [ 0.05*ones(6,1);  linspace(0.05, 2, Nframes-6)'] ;
t_delay(1:2:end) = t_delay(2:2:end);
timing_parms.t_delay        = t_delay;
timing_parms.t_adjusts      = (5 - 0.05- t_aq) * ones(Nframes,1) -t_delay;
timing_parms.labelcontrol   = isLabel;
timing_parms.labelcontrol(1:2)= -1;

timing_parms.order          = 1;
timing_parms.t_aq           = t_aq;
timing_parms.label_type     = 'FTVSI-sinc'; % 'BIR8'
timing_parms.readout_type   = 'FSE'; %'GRE';
%%%

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01;
parms.cbva =    0.01;
parms.bat =     0.1;
parms.r1tis =   1/0.9;%1.4;
parms.flip =    deg2rad(40);
parms.r2tis=    1/0.090;

doSub = 1;

parms.bat = 0.1;
s1 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);
parms.bat = 0.2;
s2 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);
parms.bat = 0.5;
s5 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);

subplot(211)
plot(s1); hold on
plot(s2)
plot(s5)
legend('100ms', '200ms', '500ms')
title('with Arterial Suppression')

%%
timing_parms.doArtSup       = 0 *ones(Nframes,1);

parms.bat = 0.1;
s1 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);
parms.bat = 0.2;
s2 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);
parms.bat = 0.5;
s5 = gen_signals_vs_220421(parms, timing_parms, 0, 1, 1e-3);

subplot(212)
plot(s1); hold on
plot(s2)
plot(s5)
title('without Arterial Suppression')
legend('100ms', '200ms', '500ms')
