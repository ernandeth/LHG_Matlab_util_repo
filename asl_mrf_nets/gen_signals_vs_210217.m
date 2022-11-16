function  obs = gen_signals_vs_210217(parms, delays, aq_parms, dofigs,doSub, dt)

% constants
if nargin<6
    dt = 1e-3;      % seconds
    %dt = 0.01;      % seconds
end
r1a = 1/1.67;

lambda = 0.9;
aqwindow = 0.650 ;
Nkz = 18;

if (nargin ==0)
    %if nargin==0
    % for testing purposes, here are some default parameters:
    f=         50 /6000;
    cbva =      0.02 ;
    bat =       0.5 ;
    r1tis =     1/1.4  ;
    
    Mtis0 =     1 ;
    flip =      90*pi/180 ; % flip angle in radians
    alpha_ai = 0.85;
    alpha_ti = 0.7;
    alpha_ts = 0.17;
    
    Nframes = 4;    

    ArtSup_delay = 0.15*ones(Nframes,1);  % delay between AS pulse and acqusition
    doSub = 0;
    dofigs = 1;
    t_tags = 0.1*ones(Nframes,1);
    t_adjusts = 2*ones(Nframes, 1);
    t_delays = 1.5*ones(Nframes,1);
    
    
    tmp =  zeros(Nframes,1);
    tmp(1:2:end) = 1;
    labelcontrol = tmp;
    labelcontrol(1:2) = -1;
    
    doArtSup = 0*ones(Nframes,1);
    

else
    t_tags = aq_parms.t_tag;  % zero if you have a single pulse
    t_tags(:) = 0;  % <---- no t_tag in the new sequence
%     t_delays = aq_parms.t_delays;
    t_delays = delays;
    t_adjusts = aq_parms.t_adjusts;
    labelcontrol = aq_parms.labelcontrol;    
    doArtSup = aq_parms.doArtSup;
    ArtSup_delay = aq_parms.ArtSup_delay ; % delay between AS pulse and acqusition
    
    aqwindow = aq_parms.t_aq(1);
    
    f = parms.f;
%     Mtis0 = parms.Mtis0;
    Mtis0 = 1;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    alpha_ai = parms.alpha_ai; % arterial inversion efficiency
    alpha_ti = parms.alpha_ti; % tissue inversion efficiency
    alpha_ts = parms.alpha_ts; % T2 weighting in tissue due to arterial suppression 
%     alpha = 0.8;
    
end
 t_tags(:) = 0;

nbat = floor(bat/dt);
Nframes = length(t_delays);
obs = zeros(Nframes,1);
%
t_tags = 0;
%
TR = t_adjusts + t_tags + t_delays + aqwindow;
begTR = cumsum(TR);
begTR = [0; begTR];  % beginning of each TR period
begTR = begTR(1:end-1);
duration = sum(TR) ;


Npts = floor((duration) / dt);
timevec = linspace(0, duration, Npts);

% indicator function for  labeling pulses 
% whether the pulse is velocity selectivity is specified by "labelcontrol"
vsfun = zeros(Npts,1);

inds = floor((begTR + t_adjusts) /dt);
inds = sort(inds);

% in the scanner file:
% labelcontrol = 1 means vs pulse
% labelcontrol = 0  means non-vs pulse
% labelcontrol = -1  means No pulse at all

% in THIS code :
% vsfun(n) = 1 means vs pulse
% vsfun(n) = -1  means non-vs pulse
% vsfun(n) = 0  means No pulse at all

tmp = zeros(size(labelcontrol));
tmp(labelcontrol==1) = 1;

vs_inds = inds .* tmp;  % vel. selective inversion
vs_inds = inds(vs_inds>0);

tmp = zeros(size(labelcontrol));
tmp(labelcontrol==0) = 1;
ns_inds = inds .* tmp;  % non-selective inversion
ns_inds = inds(ns_inds>0);

vsfun(vs_inds) = (1-alpha_ai)/2;
vsfun(ns_inds) = alpha_ai;

% indicator function for arterial suppression:
% (just add more velocity selective pulses right before acquistion)
% indicator for  arterial saturation pulses : 
asfun = zeros(Npts,1);

% indicator function for tissue T2 effects from arterial suppression:
% When we apply only the RF of arterial saturation pulses, we don't crush
% the arterial signal but we still get some T2 suppression on the tissue.
tsfun = zeros(Npts,1);


% the indices of when the AS pulses can be applied:
inds = floor((begTR + t_adjusts + t_tags + t_delays - ArtSup_delay) /dt);

% in the scanner :
% doArtSuppression = 1 means AS pulse
% doArtSuppression = -1  means No pulse at all
% doArtSuppression = 0  means non AS pulse (T2 only)

tmp = zeros(size(doArtSup));
tmp(doArtSup >=0 ) = 1;  % both zeros and ones prodice T2 effects
ts_inds = inds .* tmp;    
ts_inds = inds(ts_inds>0);
tsfun(ts_inds) = 1;

tmp = zeros(size(doArtSup));
tmp(doArtSup ==1 ) = 1;  % both zeros and ones prodice T2 effects
as_inds = inds .* tmp;    
as_inds = inds(as_inds>0);
asfun(as_inds) = 0.5;


%tmp = zeros(size(doArtSup));
%tmp(doArtSup==0) = 1;
% as_inds = inds .* tmp; 
% as_inds = inds(as_inds>0);
% asfun(as_inds) = alpha_as;

% 12/15/20: construct the arterial response  function (ARF)
tt = linspace(0,2, 2/dt);
decay = exp(-tt*r1a);
D = erf(tt/bat*2);
Arf = decay;
Arf_ex = decay.* D;

% effect of the VSI pulses on arterial magnetization:
Ma = conv(vsfun, Arf);
Ma = Ma(1:length(vsfun));

% effect of Arterial suppression pulses on arterial magnetization
Ma2 = conv(asfun, Arf);
Ma2 = Ma2(1:length(asfun));

% combine the VS and AS effects:  the  arterial suppression sets the
% magentization to zero - overrides the effect of the VSI pulses
inds = find(Ma2~=0);
Ma(inds)=Ma2(inds);

% arterial exchange compartment:  dispersion and delay 
D = exp(-(tt-bat).^2/bat*100);
Ma_ex = conv(Ma, D);  
Ma_ex = Ma_ex(1:length(Ma));
% scale to include T1 decay:
Ma_ex = Ma_ex/max(Ma_ex) * max(Ma)*exp(-bat*r1a);
%Ma = Ma/max(Ma) * max(Ma)*exp(-bat*r1a);

% Now scale the bolus concentration into magnetization:
Ma = 1-2*Ma;
Ma_ex = 1-2*Ma_ex;

% indicator for acquisition
img_sat = cos(flip)^Nkz;

img_time = floor(aqwindow / dt);  % how long it takes to collect image

aqfun = zeros(Npts,1);
inds = floor((begTR + t_adjusts + t_tags+ t_delays) /dt);
aqfun(inds) = 1;

% Magnetization time courses
M = ones(Npts, 1);
%Ma      = ones(Npts, 1) * exp((1/75-1/150)*35);%  <----  LHG 10/13/20: Need to correct for T2 weighting.
%Ma_ex   = ones(Npts, 1) * exp((1/75-1/150)*35); % 

aq = 1;
Ma0 = Ma(1);;
acqtime = -1;
batcount = 1;


%precalculate loop variables
sfp = sin(flip);

% padding added on 2020-09-25 for patching up 'array index out of bounds'
% error
aqfun = [aqfun;zeros(floor(ArtSup_delay(end)/dt),1)];
ArtSup_delay = [ArtSup_delay;ArtSup_delay(end)];
 
MM=5;
for n = nbat+MM+1:length(timevec)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The diff eqs governing the three pools of protons
    % in the voxel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  blood in artery
    %dMa = (Ma0 - Ma(n-1))*r1a ;
    %Ma(n) = Ma(n-1) + dMa*dt;
    
    % considering a flow through compartment:
    % blood in artery at exchange site separately: 
    % same thing, but has BAT lag (see below)    
    % dMa_ex = (Ma0 - Ma_ex(n-1))*r1a ;
    % Ma_ex(n) = Ma_ex(n-1) + dMa_ex*dt;
    %
    % OR ...
    % the contents of the Ma compartment feed into the exchange compartment
    % with a delay. This also means some more T1 decay before it can
    % exchange
   
    %S =  Ma(n-nbat);
    %S =  mean(Ma(n-nbat-MM:n-nbat+MM));
    
    %Ma_ex(n) = Ma0-(Ma0-S)*exp(-bat*r1a);
    
    % the modified Bloch equation has
    % t1 decay, , inflow, outflow
    dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n) - f * M(n-1)/lambda;
    M(n) = M(n-1) + dM*dt;
%     MM=5;
%     S = [M(n-MM:n-1)]' * [1:MM]'/sum(1:MM);    
%     M(n) = S + dM*dt;

    %%%%%%%%%%%%%%%%%%%
    % update the Tissue and Arterial Compartments when the pulses are
    % applied
    %%%%%%%%%%%%%%%%%%%
    
    % selective pulse
    if vsfun(n-1) ~= 0
        M(n) =  M(n) * (1-2*alpha_ti);
    end
    
    % T2 effects from the arterial selective saturation pulse
    % whether the crusher gradients are on or not.
    if tsfun(n-1) ~= 0 
        M(n) =  M(n) * abs(1-2*alpha_ts); % T2 effect only : no arterial suppression
    end
    
    
    %%%%%%%%%%%%%%%%%%%
    % Arterial Exchange Compartment:
    % pulses look like they have a lag and corresponding decay.
    % selective pulse
    %%%%%%%%%%%%%%%%%%%
    %{
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
    % LHG 2020.10.19: 
    % arterial suppression (vel non-selective) pulse - still has T2 effects
    if asfun(n-nbat) == -1
        Ma_ex(n) = Ma(n);
    end
    %}
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Readout  from the different pools
    % dependence on flip angle
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if  aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments
        % we put inMa here because it has the right relaxation due to BAT
        obs(aq) = ((1-cbva)*M(n-1) + cbva*Ma(n-1)) * sfp;
        obs_t(aq) = M(n-1) * sfp;
        obs_a(aq) = Ma(n-1) * sfp;

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
        %Ma(n) =  Ma(acqtime) * img_sat;
        %Ma_ex(n) =  Ma_ex(acqtime) * img_sat;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % influx of fresh spins:  we ASSUME that all the arterial spins are
    % fresh at the end of the acquisition.  
    % Reset the matnetization after img_time.
    % Make sure t_adjust > 1.8 for this to hold in reality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (acqtime > 0) && ((n-acqtime) == img_time)
        %Ma(n) =  Ma0;
        %Ma_ex(n) =  1;
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
    
    %some times they don't match up.  flooring error somewhere?
    ed = length(M);
    labelfun = vsfun;
    labelfun = labelfun(1:length(timevec));
    aqfun = aqfun(1:length(timevec));
    
    %figure()
    subplot(311)
    stem(timevec(1:ed),labelfun(1:ed),'r')
    hold on
    stem(timevec(1:ed),aqfun(1:ed),'k');
    stem(timevec(1:ed),asfun(1:ed),'m');
    
    legend('Label','Aq', 'ArtSup')
    hold off
    axis tight
    
    subplot(312)
    plot(timevec(1:ed),M(1:ed))
    hold on
    plot(timevec(1:ed),Ma(1:ed),'r')
    plot(timevec(1:ed),Ma_ex(1:ed),'g')
    plot(timevec(1:ed),aqfun(1:ed),'k--')
    legend('Tis', 'Art', 'ArtX')
    
    hold off
    axis tight
    grid on
    
    subplot(313)
    plot(obs)
    hold on
    plot(obs_t)
    %plot(cbva*obs_a(2:end))
    legend('Both', 'Tis') % , 'Art')
    hold off
    grid on
    
    drawnow()
end



return
