% creating ASL timing schedule for CBF measurement

Nframes = 16;
t_aq = 0.57
t_aq = 0.5

isLabel     = ones(Nframes,1);
isLabel(1:2:end) = 0;

timing_parms.doArtSup       = ones(Nframes,1);
timing_parms.doArtSup(1:2)  = 0;
timing_parms.ArtSup_delay   = 0.05 * ones(Nframes,1);
timing_parms.Nlabel_group   = 1;
timing_parms.t_tag          = 0.3 *ones(Nframes,1 );
timing_parms.t_delay        = 1.3 * ones(Nframes,1) ;
timing_parms.t_adjusts       = (5-1.3 - 0.05- t_aq) * ones(Nframes,1) ;
timing_parms.labelcontrol   = isLabel;
timing_parms.labelcontrol(1:2)= -1;

timing_parms.order          = 1;
timing_parms.t_aq           = t_aq;
timing_parms.label_type     = 'BIR8'; % 'BIR8'
timing_parms.readout_type   = 'FSE'; %'GRE';
        
make_schedule_files(timing_parms, 'cbf_vss_schedule')

%% creating ASL timing schedule for BAT measurement
t_aq        = 0.51;

Nframes = 50;


isLabel     = ones(Nframes,1);
isLabel(1:2:end) = 0;

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
timing_parms.label_type     = 'BIR8'; % 'BIR8'
timing_parms.readout_type   = 'FSE'; %'GRE';
        
make_schedule_files(timing_parms, 'bat_vss_schedule_phantom_noAS')

quick_mrf_sensitivity(timing_parms)

%% creating ASL timing schedule for multi-PLD measurement

Nframes = 32;

t_aq        = 0.51;
isLabel     = ones(Nframes,1);
isLabel(1:2:end) = 0;

timing_parms.doArtSup       = ones(Nframes,1);
timing_parms.doArtSup(:)  = 0;

timing_parms.ArtSup_delay   = 0.05 * ones(Nframes,1);
timing_parms.Nlabel_group   = 1;
timing_parms.t_tag          = 0.3 *ones(Nframes,1 );

t_delay  = linspace(0.1, 3, Nframes)' ;
t_delay(1:2:end) = t_delay(2:2:end);

timing_parms.t_delay        = t_delay;
timing_parms.t_adjusts      = (5 - 0.05- t_aq) * ones(Nframes,1) -t_delay;
timing_parms.labelcontrol   = isLabel;
%timing_parms.labelcontrol(1:2)= -1;

timing_parms.order          = 1;
timing_parms.t_aq           = t_aq;
timing_parms.label_type     = 'BIR8'; % 'BIR8'
timing_parms.readout_type   = 'FSE'; %'GRE';
        
make_schedule_files(timing_parms, 'multiPLD_vss_schedule')


quick_mrf_sensitivity(timing_parms)

%% creating ASL timing schedule for CBV measurement

Nframes = 20;


isLabel     = ones(Nframes,1);
isLabel(1:2:end) = 0;

timing_parms.doArtSup       = ones(Nframes,1);
timing_parms.doArtSup(1:2:end)  = 0;
timing_parms.doArtSup(1:2)  = 0;

timing_parms.ArtSup_delay   = 0.05 * ones(Nframes,1);
timing_parms.Nlabel_group   = 1;
timing_parms.t_tag          = 0.3 *ones(Nframes,1 );
timing_parms.t_delay        = 1.3 * ones(Nframes,1) ;
timing_parms.t_adjusts      = (5-1.3 - 0.05- t_aq) * ones(Nframes,1) ;
timing_parms.labelcontrol   = isLabel;
timing_parms.labelcontrol(:)= -1;

timing_parms.order          = 1;
timing_parms.t_aq           = t_aq;
timing_parms.label_type     = 'BIR8'; % 'BIR8'
timing_parms.readout_type   = 'FSE'; %'GRE';
        
make_schedule_files(timing_parms, 'cbv_vss_schedule')
quick_mrf_sensitivity(timing_parms)
