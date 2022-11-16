clear all
%
%addpath  /home/hernan/matlab/flow/fingerprint
read_parms_file = 0
create_schedule_files = 1;
schedule_name = 'vsi_lindelay_100_noartsup'
    
dt = 1e-4; % <----  note 


% create ASL sequence timing
if read_parms_file
    % read in schedule from file:
    aq_parms.t_tag          = load('t_tags.txt');
    aq_parms.t_delays       = load('t_delays.txt');
    aq_parms.t_adjusts      = load('t_adjusts.txt');;
    aq_parms.labelcontrol   = load('isVelocitySelective.txt');
    aq_parms.doArtSup       = load('doArtSuppression.txt');;
    aq_parms.ArtSup_delay   = 0.1 *ones(size(delays)); % delay between AS pulse and acqusition
    aq_parms.t_aq           = load('t_aqs.txt');
    Nframes = length(t_delay);
else
    % incrementing schedule:
    %
    Nframes = 100;
    idx = randperm(Nframes);
    TR = 4.5;
    aq_parms.labelcontrol           = ones(Nframes,1);
    aq_parms.labelcontrol(2:2:end)  = 0;
    aq_parms.labelcontrol(1:4)      = -1;
     
    
    aq_parms.doArtSup               = ones(Nframes,1);
     aq_parms.doArtSup(:)            = 0;
     aq_parms.doArtSup(:)            = -1;
    
    aq_parms.t_aq                   = 0.65 *ones(Nframes,1);
    aq_parms.ArtSup_delay           = 0.1 *ones(Nframes,1); % delay between AS pulse and acqusition
    aq_parms.t_adjusts              = 1.5 *ones(Nframes,1);
    aq_parms.t_tag                  = zeros(Nframes,1);
    aq_parms.t_delays               =  linspace( 0.1, 1.5, Nframes)';
    aq_parms.t_delays(2:2:end)      =  aq_parms.t_delays(1:2:end);
    aq_parms.t_adjusts              = TR - aq_parms.t_delays;


end

% default parameters
parms.f         = 60/6000;
parms.cbva      = 0.01 ;%0.02;
parms.bat       = 0.5;
parms.r1tis     = 1/1.4;

% fixed defaults
parms.flip      = deg2rad(90);
parms.alpha_ai  = 0.8; % arterial inversion efficiency
parms.alpha_ti  = 0.75; % tissue inversion efficiency
parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression
parms.Mtis0     = 1;

%%
% sanity checks: show the signal for one case
doSub = 1;
doFigs = 0;

parms.f = 60/6000;
parms.bat = 0.8;
parms.cbva = 0.01;
entry1 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
parms.bat = 0.4;
entry2 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));

figure
subplot(211)
plot(entry1); 
hold on
plot(entry2); 
title('subtraction signals')
legend('BAT = 0.8', 'BAT = 0.4')

subplot(212)
plot(entry2-entry1)
title('Difference')

%%
parms.f = 60/6000;
parms.bat = 0.4;
parms.cbva = 0.01;
entry1 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
parms.f = 30/6000;
entry2 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));

figure
subplot(211)
plot(entry1); 
hold on
plot(entry2); 
title('subtraction signals')
legend('CBF = 60', 'CBF = 30')

subplot(212)
plot(entry2-entry1)
title('Difference')
%%

parms.f = 60/6000;
parms.bat = 0.4;
parms.cbva = 0.01;
entry1 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
parms.cbva = 0.02;
entry2 = (gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));

figure
subplot(211)
plot(entry1); 
hold on
plot(entry2); 
title('subtraction signals')
legend('CBV = 0.01', 'CBV = 0.02')

subplot(212)
plot(entry2-entry1)
title('Difference')
%%


if create_schedule_files
    
    mkdir (schedule_name)
    cd(schedule_name)
    
    tmp = aq_parms.t_tag;
    save t_tags.txt tmp -ascii
    
    tmp = aq_parms.t_delays;
    save t_delays.txt tmp -ascii
    
    tmp =aq_parms.t_adjusts;
    save t_adjusts.txt tmp -ascii
    
    tmp = aq_parms.labelcontrol;
    save isVelocitySelective.txt tmp -ascii
    
    tmp= aq_parms.doArtSup ;       
    save doArtSuppression.txt tmp -ascii
   
    tmp = aq_parms.t_aq;
    save t_aqs.txt tmp -ascii;
    
    cd ..
end