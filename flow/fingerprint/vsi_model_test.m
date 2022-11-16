clear all
%
addpath  /home/hernan/matlab/flow/fingerprint
read_parms_file = 0

dt = 1e-3; % <----  note 


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
    % alternative schedule:
    %
    Nframes = 100;
    idx = randperm(Nframes);
    
    aq_parms.labelcontrol           = ones(Nframes,1);
    aq_parms.labelcontrol(2:2:end)  = 0;
    % aq_parms.labelcontrol(idx(1:end/2))  = 0;    % randomize labels
    % aq_parms.labelcontrol(1:8)        = -1;
     
    
    aq_parms.doArtSup               = ones(Nframes,1);
    % aq_parms.doArtSup(5:5:end/2)    = 0;
    % aq_parms.doArtSup(1:4)          = -1;
    
    aq_parms.t_aq                   = 0.5 *ones(Nframes,1);
    aq_parms.ArtSup_delay           = 0.1 *ones(Nframes,1); % delay between AS pulse and acqusition
    aq_parms.t_adjusts              = 1.5 *ones(Nframes,1);
    aq_parms.t_tag                  = zeros(Nframes,1);
    %aq_parms.t_delays               =  0.6 + 0.9*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    aq_parms.t_delays               =  1.2 + 0.3*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    aq_parms.t_delays               =  0.1 + 1.1*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    %aq_parms.t_delays(2:2:end)      = aq_parms.t_delays(1:2:end);


end

% default parameters
parms.f         = 60/6000;
parms.cbva      = 0.02;
parms.bat       = 0.5;
parms.r1tis     = 1/1.4;

% fixed defaults
parms.flip      = deg2rad(90);
parms.alpha_ai  = 0.8; % arterial inversion efficiency
parms.alpha_ti  = 0.75; % tissue inversion efficiency
parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression
parms.Mtis0     = 1;

%
% sanity check: show the signal for one case
doSub = 0;
doFigs = 0;
figure
subplot(211)
parms.f = 60/6000;
parms.bat = 0.2;
entry1 = abs(gen_signals_vs_201020(parms, aq_parms, doFigs,doSub, dt));
entry2 = abs(gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
plot(entry1/entry1(1)); 
hold on
plot(entry2/entry2(1)); 

subplot(212)
plot(entry1-entry2)
