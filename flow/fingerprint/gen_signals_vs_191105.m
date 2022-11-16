function  obs = gen_signals_vs_191105(parms, aq_parms, dofigs,doSub)

% constants
dt = 0.0001;      % seconds
% dt = 0.01;      % seconds
r1a = 1/1.67;

lambda = 0.9;
aqwindow = 0.650 ;
Nkz = 16;

if (nargin ==0)
    %if nargin==0
    % for testing purposes, here are some default parameters:
    f=         50 /6000;
    Mtis0 =     1 ;
    cbva =      0.02 ;
    bat =       0.2 ;
    r1tis =     1/1.4  ;
    flip =      90*pi/180 ; % flip angle in radians
    
    alpha_ai =     0.81;  % efficiency of VSI pulse on the blood
    alpha_ti =     0.66;  % efficiency of the VSI pulse on the tissue
    alpha_ts =     0.77;  % efficiency of the VSSat pule on the tissue
    
    Nframes = 4;    

    ArtSup_delay = 0.5;  % delay between AS pulse and acqusition
    doSub = 1;
    dofigs = 1;
    t_tags = 1.5*ones(Nframes,1);
    t_adjusts = 1*ones(Nframes, 1);
    t_delays = linspace(0.15, 3, Nframes)';
    t_delays(1:2:end)  = t_delays(2:2:end);
    t_delays = 1.8*ones(Nframes,1);
    
    tmp =  zeros(2*Nframes,1);
    tmp(3:4:end) = 1;
    tmp(4:4:end) = 1;
    labelcontrol = tmp;
    
    tmp =  zeros(Nframes,1);
    tmp(1:2:end) = 1;
    labelcontrol = tmp;
    
    doArtSup = ones(Nframes,1);
    doArtSup(:) = 1;

else
    t_tags = aq_parms.t_tags;  % zero if you have a single pulse
    t_delays = aq_parms.t_delays;
    t_adjusts = aq_parms.t_adjusts;
    labelcontrol = aq_parms.labelcontrol;    
    doArtSup = aq_parms.doArtSup;
    ArtSup_delay = aq_parms.ArtSup_delay ;% delay between AS pulse and acqusition
    alpha = aq_parms.alpha;
    
    f = parms.f;
    Mtis0 = parms.Mtis0;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    
end


nbat = round(bat/dt);
Nframes = length(t_delays);
obs = zeros(Nframes,1);

TR = t_adjusts + t_tags + t_delays + aqwindow;
begTR = cumsum(TR);
begTR = [0; begTR];  % beginning of each TR period
begTR = begTR(1:end-1);
duration = sum(TR) ;


Npts = round((duration) / dt);
timevec = linspace(0, duration, Npts);

% indicator for  labeling pulses : a pair of pulses separated by t_tag
% whether the pulse is velocity selectivity is specified by "labelcontrol"
vsfun = zeros(Npts,1);

% indicator for  arterial saturation pulses : 
asfun = zeros(Npts,1);


inds = round((begTR + t_adjusts) /dt);
%inds2 = round((begTR + t_adjusts + t_tags) /dt);% (second pulse here)
%inds = [inds; inds2];
inds = sort(inds);

vs_inds = inds .* labelcontrol;  % vel. selective
vs_inds = inds(inds>0);
vsfun(vs_inds) = 1;

ns_inds = inds .* (~labelcontrol);  % non-selective
ns_inds = inds(ns_inds>0);
vsfun(ns_inds) = -1;

% indicator for arterial suppression:
% (just add more velocity selective pulses right before acquistion)
inds = round((begTR + t_adjusts + t_tags+ t_delays - ArtSup_delay) /dt);
avs_inds = inds .* doArtSup;  % vel. selective
avs_inds = inds(avs_inds>0);
asfun(avs_inds) = 1;


% indicator for acquisition
img_sat = cos(flip)^Nkz;
img_time = round(aqwindow / dt);  % how long it takes to collect image

aqfun = zeros(Npts,1);
inds = round((begTR + t_adjusts + t_tags+ t_delays) /dt);
aqfun(inds) = 1;

% Magnetization time courses
M = ones(Npts, 1);
Ma = ones(Npts, 1);
Ma_ex = ones(Npts, 1);

aq = 1;
Ma0 = 1;
acqtime = -1;
batcount = 1;

for n = nbat+1:length(timevec)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The diff eqs governing the three pools of protons
    % in the voxel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  blood in artery
    dMa = (Ma0 - Ma(n-1))*r1a ;
    Ma(n) = Ma(n-1) + dMa*dt;
    
    % blood in artery at exchange site separately: 
    % same thing, but has BAT lag (see below)    
    dMa_ex = (Ma0 - Ma_ex(n-1))*r1a ;
    Ma_ex(n) = Ma_ex(n-1) + dMa_ex*dt;
    
    
    % the modified Bloch equation has
    % t1 decay, , inflow, outflow
    dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n)  - f * M(n-1)/lambda;
    M(n) = M(n-1) + dM*dt;
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Readout  from the different pools
    % dependence on flip angle
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if  aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments
        % we put inMa here because it has the right relaxation due to BAT
        obs(aq) = ((1-cbva)*M(n-1) + cbva*Ma(n-1)) * sin(flip);
        obs_t(aq) = M(n-1) * sin(flip);
        obs_a(aq) = Ma(n-1) * sin(flip);

        aq = aq +1;
        acqtime = n;  % mark the begining of the acquisition
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % readout causes partial saturation of spins
    % readout (acqusition sequence)
    % flip angle dependence calculated in img_sat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (acqtime > 0) && ((n-acqtime) < img_time)
        M(n ) =   M(acqtime) * img_sat;
        Ma(n) =  Ma(acqtime) * img_sat;
        Ma_ex(n) =  Ma_ex(acqtime) * img_sat;
    end
    
    %%%%%%%%%%%%%%%%%%%
    % update the Tissue and Arterial Compartments when the pulses are
    % applied
    %%%%%%%%%%%%%%%%%%%
    
    % selective pulse
    if vsfun(n-1) == 1
        Ma(n) = Ma(n)*(1-2*alpha_ai);
        M(n) =  M(n) * alpha_ti;
    end
    
    % non-selective pulse
    if vsfun(n-1) == -1
        Ma(n) = Ma(n) * (1-2*alpha_ai);
        M(n)  = M(n) *  (1-2*alpha_ti);
    end
    
     % arterial selective saturation pulse
    if asfun(n-1) == 1
        Ma(n) = 0;
        M(n) =  M(n) * alpha_ts;
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Arterial Exchange Compartment:   
    % pulses look like they have a lag and corresponding decay.
    % selective pulse
    %%%%%%%%%%%%%%%%%%%
    if n >= nbat
        if vsfun(n-nbat) == 1
            Ma_ex(n) = Ma(n);
        end
        % non-selective pulse
        if vsfun(n-nbat) == -1
            Ma_ex(n) = Ma(n);
        end
        % arterial suppression (vel selective) pulse
        if asfun(n-nbat) == 1
            Ma_ex(n) = Ma(n);
        end
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % influx of fresh spins:  we ASSUME that all the arterial spins are
    % fresh at the end of the acquisition.  
    % Reset the matnetization after img_time.
    % Make sure t_adjust > 1.8 for this to hold in reality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (acqtime > 0) && ((n-acqtime) == img_time)
        Ma(n) =  1;
        Ma_ex(n) =  1;
    end
    
end

if doSub
    obs_a = (obs_a(1:2:end) - obs_a(2:2:end));
    obs_t = (obs_t(1:2:end) - obs_t(2:2:end));
    obs =   (obs(1:2:end) -   obs(2:2:end)); %./ ...
       % (obs(1:2:end) +   obs(2:2:end)/2);
end

obs = obs';

if dofigs
    
    %some times they don't match up.  rounding error somewhere?
    
    labelfun = vsfun;
    labelfun = labelfun(1:length(timevec));
    aqfun = aqfun(1:length(timevec));
    
    %figure()
    subplot(311)
    stem(timevec,labelfun,'r')
    hold on
    stem(timevec,aqfun,'g');
    stem(timevec,0.5*asfun,'m');
    
    legend('Label','Aq', 'ArtSup')
    hold off
    axis tight
    
    subplot(312)
    plot(timevec,M)
    hold on
    plot(timevec,Ma,'r')
    plot(timevec,Ma_ex,'m')
    plot(timevec,aqfun,'g--')
    legend('Tis', 'Art', 'ArtX')
    
    hold off
    axis tight
    grid on
    
    subplot(313)
    plot(obs(2:end))
    hold on
    plot(obs_t(2:end))
    %plot(cbva*obs_a(2:end))
    legend('Both', 'Tis') % , 'Art')
    hold off
    grid on
end



return
