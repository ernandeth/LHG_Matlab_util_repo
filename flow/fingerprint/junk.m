figure

phys_parms.f =         0.01;
phys_parms.mtis0 =     1 ;
phys_parms.cbva =      0.03 ;
phys_parms.transit=    1 ;
phys_parms.kfor =     0.2 ;
phys_parms.r1tis =     1/1.4  ;
phys_parms.beta =      90*pi/180 ; % flip angle in radians
phys_parms.L = 1;
phys_parms.Disp =      35;
phys_parms.Ptime =     0.5;


psdseqtime = 0.0346;
nslices = 1;
%
for p=[-0.002:0.001:0.002]
    phys_parms.f =      0.05 + p ;

    label_dur = load('t_tags.txt');
    PID = load('t_delays.txt');
    % PID = PID + 0.014; % account for the fatsat pulse
    
    % PID = PID + p;   % debugging ... ?
    
    TR = label_dur + PID + nslices*psdseqtime + 0.035;
    
   % TR = label_dur + PID + nslices*psdseqtime ;  % no delay after acq.
    
    
    %% testing only :
    % label_dur = sort(repmat(linspace(0.1,3,50),1,2))';
    % PID =1.2 * ones(size(label_dur));
    %
    %  label_dur= 3.2* ones(6,1);
    %  PID = 1.8 * ones(6,1);
    %  TR = 5.5 * ones(6,1);
    %
    %%
    aqtimes = label_dur + PID ;  % acquisition happens 10 ms after labeling
    duration = sum(TR)+2;
    
    timing_parms.PID = PID;
    timing_parms.label_dur = label_dur;
    timing_parms.TR = TR;
    
    % Just a quick test to make sure everything works as expected
    obs = gen_signals_140420(phys_parms, timing_parms, 0,0);
    obs = obs - mean(obs);
    obs = obs / norm(obs);
    subplot(211)
    plot(obs); hold on
    
    subplot(212)
    plot(obs(1:2:end) - obs(2:2:end)); hold on
    
    drawnow
    
end