function  obs = gen_signals_vs_211024(parms, delays, aq_parms, dofigs,doSub, dt)

% constants
if nargin<6
    dt = 1e-3;      % seconds
    dt = 0.01;      % seconds
end
r1a = 1/1.67;

lambda = 0.9;
aqwindow = 0.650 ;
Nkz = 18;

r2art  = 1/0.150 ;
eTE_ftvsi = 0.0313;     % effective echo time of FT-VSI sinc-pulse
eTE = 0.0207;           % effective echo time of sBIR8-VSS pulse

if (nargin ==0)
    %if nargin==0
    % for testing purposes, here are some default parameters:
    f=         100 /6000;
    cbva =      0.02 ;
    bat =       0.2 ;
    r1tis =     1/1.4  ;
    Mtis0 =     1 ;
    flip =      40*pi/180 ; % flip angle in radians
    
    r2tis = 1/0.090 ;
        
    Nframes = 4;
    
    ArtSup_delay = 0.2*ones(Nframes,1);  % delay between AS pulse and acqusition
    doSub = 0;
    dofigs = 1;
    t_tags = 0.1*ones(Nframes,1);
    t_adjusts = 2.5*ones(Nframes, 1);
    t_delays = 1.3*ones(Nframes,1);
    
    t_adjusts = 1.5*ones(Nframes, 1);
    t_delays = 0.5*ones(Nframes,1);
    
    tmp =  zeros(Nframes,1);
    %tmp(2:2:end) = 1;
    labelcontrol = tmp;
    %labelcontrol(1) = -1;
    order = 1;
    
    doArtSup =  ones(Nframes,1);
    doArtSup(2:2:end) = 0;

    label_type = 'FTVSI-sinc'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    label_type = 'BIR8';
    %label_type = 'BIR8inv';
    readout_type = 'GRE'; 
    readout_type = 'FSE'; 

else
    t_tags = aq_parms.t_tag;  % zero if you have a single pulse
    t_tags(:) = 0;  % <---- no t_tag in the new sequence
    %   t_delays = aq_parms.t_delays;
    t_delays = delays;
    t_adjusts = aq_parms.t_adjusts;
    labelcontrol = aq_parms.labelcontrol;
    doArtSup = aq_parms.doArtSup;
    ArtSup_delay = aq_parms.ArtSup_delay ; % delay between AS pulse and acqusition
    
    aqwindow = aq_parms.t_aq(1);
    label_type = aq_parms.label_type;
    order = aq_parms.order;
    readout_type = aq_parms.readout_type;
    
    
    f = parms.f;
    %   Mtis0 = parms.Mtis0;
    Mtis0 = 1;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    r2tis = parms.r2tis;
    
    
end

t_tags = round(100*t_tags)/100;
t_delays = round(100*t_delays)/100;
t_adjusts = round(100*t_adjusts)/100;
ArtSup_delay = round(100*ArtSup_delay)/100;
aqwindow = round(100*aqwindow)/100;

% LHG 3/1/2022
% in the pulse sequence, the artsup_delay happens after the t_delay
% but that's not the case here.  We must account for it
t_delays = t_delays + ArtSup_delay;

%------
% calculate effects of VS pulses given tissue T2
% alpha is the degree of inversion imparted by the pulse
T2loss     = exp(-eTE*r2tis);
T2loss_art = exp(-eTE*r2art);

T2loss_ftvsi     = exp(-eTE_ftvsi*r2tis);
T2loss_art_ftvsi = exp(-eTE_ftvsi*r2art);

switch(label_type)
    case 'BIR8inv'
        % if the first pulse is  BIR8-inv : inversion efficiency with  T2 effects
        % tissue gets inverted too:  M(n) =  M(n-1) * (-alpha1_tis_sel);
        %
        alpha1_tis_sel = -T2loss;
        alpha1_tis_ns  = -T2loss;
        alpha1_art_sel = 0.5;
        alpha1_art_ns  = (1+T2loss_art)/2;
        %}
        
    case 'BIR8'
        % first pulse BIR8-sat : T2 effects
        % tissue does this: M(n) =  M(n-1) * alpha1_tis_sel
        %
        alpha1_tis_sel = T2loss;
        alpha1_tis_ns  = T2loss;
        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
        alpha1_art_sel = 0.5;
        alpha1_art_ns  = (1-T2loss_art)/2;
        %
    case 'FTVSI-sinc'
        % first pulse BIR8-sat : T2 effects
        % tissue does this: M(n) =  M(n-1) * alpha1_tis_sel
        %
        alpha1_tis_sel = -T2loss_ftvsi;
        alpha1_tis_ns  = -T2loss_ftvsi;
        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
        alpha1_art_sel = (1+T2loss_art_ftvsi)/2;
        alpha1_art_ns  = (1-T2loss_art_ftvsi)/2;
end

% second pulse : T2 effects only on the tissue
alpha2_tis_sel = T2loss;
alpha2_tis_ns  = T2loss;
% second pulse : SATURATion or T2 effects only on the artery
alpha2_art_sel = 0.5;
alpha2_art_ns  = (1-T2loss_art)/2;
%------
tmplabelcontol= labelcontrol;
if order==2
    tmplabelcontrol(labelcontrol==0) = 1;
    tmplabelcontrol(labelcontrol==1) = 0;
    labelcontrol = tmplabelcontrol;
end

%---
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

vs_inds = inds(:) .* tmp(:);  % vel. selective inversion
vs_inds = inds(vs_inds>0);

tmp = zeros(size(labelcontrol));
tmp(labelcontrol==0) = 1;
ns_inds = inds(:) .* tmp(:);  % non-selective inversion
ns_inds = inds(ns_inds>0);

vsfun(vs_inds) = -1;% (1-alpha_art_ns)/2;
vsfun(ns_inds) = 1; % alpha_art_ns;

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
tmp(doArtSup ==1 ) = 1;  % both zeros and ones produce T2 effects
as_inds = inds .* tmp;
as_inds = inds(as_inds>0);
asfun(as_inds) = 1;

tmp = zeros(size(doArtSup));
tmp(doArtSup ==0 ) = 1;  % both zeros and ones produce T2 effects
as_inds = inds .* tmp;
as_inds = inds(as_inds>0);
asfun(as_inds) = -1;

% make an indicator for the beginning of readout
aqfun = zeros(Npts,1);
inds = floor((begTR + t_adjusts + t_tags+ t_delays) /dt);
aqfun(inds) = 1;

%------------
% LHG :2.5.22 - adapt the effect of readout for both FSE and GRE
% approximate the  effects of the read out on the Mz of tissue
% what is the "steady state" of the magnetization after the readout pulses.
if readout_type == 'FSE'
    flip = deg2rad(120);
end
cosflip = cos(flip);  % pre-compute this once to save time
%-------------

img_time = floor(aqwindow / dt);  % how long it takes to collect image
slice_time = floor(img_time/Nkz);


MM=5;

% Create arterial input function:
% Note: at the end of each bolus, we restore the magentization with fresh spins
% ie -  skip calculating the decay in the loop for that
% time point and reset the label to zero instead

bolus_length = 2; % bolus from labeling
bolus_length_aq = floor(aqwindow*100)/100; %0.5; % bolus created by aquisition
%(assume a labeling pulse is always applied after readout, to reset the arterial magnetication)
BN = floor(bolus_length/dt);
BN_aq = floor(bolus_length_aq/dt);

tt = linspace(0,bolus_length, BN);
tt_aq = linspace(0,bolus_length_aq, BN_aq);  % bolus duration from imaging

% assuming fresh blood :
bolus = exp(-tt*r1a)';
bolus_aq = exp(-tt_aq*r1a)';
Ma = zeros(Npts, 1);

for n = MM:length(timevec)
    
    if vsfun(n-1) == 1  % gradients ON
        Ma(n+1 : n+BN)= alpha1_art_sel * bolus;
    end
    
    if vsfun(n-1) == -1 % gradients OFF
        Ma(n+1:n+BN) =  alpha1_art_ns * bolus;
    end
    
    % ArtSup pulses cause saturation of the blood.
    % saturation recovery curve
    if asfun(n-1) == 1
        Ma(n+1 : n+BN)= alpha2_art_sel * bolus;
    end
    
    % Assume that Artsup control pulses don't do anything to  arterial
    % blood.  T2 effects on tissue spins only.
    % ArtSup VS control pulses cause T2 attenuation , then T1 decay
    if asfun(n-1) == -1
        if Ma(n-1) == 0
            Ma(n+1 : n+BN)= alpha2_art_ns * bolus;
        end
    end
    
    % LHG: 2/6/2022  adapt for FSE vs. GRE 
    % Finally effect of the readout is that it destroys magnetization by the end
    % acts as a shorter volume of label
    %
    if readout_type == 'FSE'
        if n>img_time+1
            if aqfun(n)==1
                %Ma(n+1 : n + BN_aq)= 0.5 * bolus_aq;
                Ma(n+1  : n+BN_aq)= 0.5 ;
                %Ma(n+img_time +1 : n+img_time+BN_aq)= 0.5 * bolus_aq;
            end
        end
    end
    %}
end



%Ma = Ma2;
% arterial exchange compartment:  dispersion and delay
%D = exp(-(tt-bat).^2/bat*100);

% change to gamma function (1/31/22)
D = tt.* exp(-(tt/bat*2).^2);
D = D / sum(D);


Ma_ex = conv(Ma, D);
Ma_ex = Ma_ex(1:length(Ma));

% scale to include T1 decay:
Ma_ex = Ma_ex/max(Ma_ex) * max(Ma)*exp(-bat*r1a);
% convert from to concentration to Magnetization
Ma_ex = 1-2*Ma_ex;
Ma = 1-2*Ma;



% Tissue Magnetization time courses
M = ones(Npts, 1);

aq = 1;
Ma0 = Ma(1);
acqtime = 0;
kz = 0;
batcount = 1;


%precalculate loop variables
sinflip = sin(flip);

% padding added on 2020-09-25 for patching up 'array index out of bounds'
% error
aqfun = [aqfun; zeros(floor(ArtSup_delay(end)/dt),1)];
ArtSup_delay = [ArtSup_delay;ArtSup_delay(end)];


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
    if vsfun(n-1) == 1
        M(n) =  M(n) * alpha1_tis_sel;
    end
    % non-selective pulse
    if vsfun(n-1) == -1
        M(n) =  M(n) * alpha1_tis_ns;
    end
    
    % T2 effects from the arterial selective saturation pulse
    % whether the crusher gradients are on or not.
    if asfun(n-1) ~= 0
        M(n) =  M(n) * (alpha2_tis_ns); % T2 effect only : no inversion
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Readout  from the different pools
    % dependence on flip angle
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if  aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments
        
        %{
        obs(aq) = ((1-cbva)*M(n) + cbva*Ma(n)) * sinflip;
        obs_t(aq) = M(n) * sinflip;
        obs_a(aq) = Ma(n) * sinflip;
        %}
        %
        obs(aq) = ((1-cbva)*M(n-1) + cbva*Ma(n-1)) * sinflip;
        obs_t(aq) = M(n-1) * sinflip;
        obs_a(aq) = Ma(n-1) * sinflip;
        %}
        
        aq = aq +1;
        acqtime = n;  % mark the begining of the acquisition
               
        M(n) = M(n) * cosflip;

        if readout_type=='FSE'
            M(n) = 0;
        end
        kz = 1;

    end


    % LHG:2/6/22 - adapted for FSE and GRE cases 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % readout causes partial saturation of spins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if  (n == acqtime + slice_time*kz)
       M(n) = M(n) * cosflip;
       kz = kz + 1;
       if kz > Nkz
           kz = 0;
           acqtime = 0;
       end
    end
    %}
    
    
   
    
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
    plot(timevec(1:ed),M(1:ed),'b')
    hold on
    plot(timevec(1:ed),Ma(1:ed),'g')
    plot(timevec(1:ed),Ma_ex(1:ed),'r')
    plot(timevec(1:ed),aqfun(1:ed),'k--')
    legend('Tis', 'Art', 'ArtX')
    axis([min(timevec) max(timevec) -1 1])
    
    hold off
  %  axis tight
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
