function  obs = gen_signals_vs_230918(parms, aq_parms, dofigs,doSub, dt)
% this version is modified to agree with the timings of the new pulse
% sequence asl3dflex
% it also considers the effect of multiple pulses on the arterial
% compartment - splits into multiple populations

% constants
if nargin<5
    dt = 1e-3;      % seconds
    %dt = 5e-3;
    %dt = 0.01;      % seconds
end

r1a = 1/1.67;
lambda = 0.9;
r2art  = 1/0.150 ;

eTE_ftvsi = 0.0313;     % effective echo time of FT-VSI sinc-pulse
eTE = 0.0207;           % effective echo time of sBIR8-VSS pulse

BGS0_length = round(2.5e-3 /dt); % (s)  duration of the global saturation pulse - magnetization reset.

if (nargin ==0)
    % for testing purposes, here are some default tissue parameters:
    % Tissue parms:
    f=          0.01 ;
    cbva =      0.02;
    bat =       0.2 ;
    r1tis =     1/1.4  ;
    Mtis0 =     1 ;
    flip =      10*pi/180 ; % flip angle in radians
    r2tis =     1/0.090 ;
    
    % include a B1 error term in the labeling pulses
    % 0 error means multiplying by 1.
    % 1% error mean multiplying by 1.01
    b1err =     0.05;

    % processing defaults
    doSub = 1;
    dofigs = 1;
        
    % pulse sequence parameters
    Nframes = 4;
    label_type =    'FTVSI-sinc'; %'FTVSS'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    RO_type =       'FSE'; %'FSE';   % 'GRE'
    t_tags =        0;% 0.1*ones(Nframes,1);
    del1 =          2*ones(Nframes, 1);
    del2 =          1.1*ones(Nframes,1);
    del3 =          0.5*ones(Nframes,1);  % delay between AS pulse and acqusition
    labelcontrol =  zeros(Nframes,1);
    labelcontrol(2:2:end)= 1;
    %labelcontrol(1:2) =   -1;
    order =         1;    
    doArtSup =      ones(Nframes,1);
    %doArtSup(end/2+1:end)= 0;
    doArtSup(1:2)=    0;
    % test case
    %{
    del1 = 1.5*ones(Nframes, 1);
    del2 = 0.5*ones(Nframes,1);
    %}
    
    RO_time = 0.750 ;  % duration of the whole readout
    Nkz = 16;
    
else
    Nkz = 16;

    % pulse sequence parms:
    t_tags =        aq_parms.t_tags;  % zero if you have a single pulse
    t_tags(:) =     0;  % <---- no t_tag in the new sequence
    del1 =          aq_parms.del1;
    del2 =          aq_parms.del2;
    del3 =          aq_parms.del3 ; % delay between AS pulse and acqusition
    labelcontrol =  aq_parms.labelcontrol;
    label_type =    aq_parms.label_type;
    order =         aq_parms.order;
    doArtSup =      aq_parms.doArtSup;
    RO_time =       aq_parms.t_aq(1);
    RO_type =       aq_parms.RO_type;
    
    % Tissue parms
    f = parms.f;
    Mtis0 = parms.Mtis0;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    flip =  parms.flip;
    r2tis = parms.r2tis;

    % include a B1 error term in the labeling pulses
    % 0 error means multiplying by 1.
    % 1% error mean multiplying by 1.01
    b1err = parms.b1err;

end

% Translate labeling error (B1err) into a scaling factor
b1err= b1err +1;

% check the order.  The 'order' flag is a quick way to reverse the order
% ie - VS becomes non-VS and non-VS becomes VS
tmplabelcontol= labelcontrol;
if order==2  %  2023.08.08....original should be  2 (just checking to see if there is an error)
    tmplabelcontrol(labelcontrol==0) = 1;
    tmplabelcontrol(labelcontrol==1) = 0;
    tmplabelcontrol(labelcontrol==-1) = -1;
    
    labelcontrol = tmplabelcontrol;
end

% rounding units
t_tags = round(100*t_tags)/100;
del1 = round(100*del1)/100;
del2 = round(100*del2)/100;
del3 = round(100*del3)/100;
RO_time = round(100*RO_time)/100;

%------
% calculate effects of VS pulses given tissue T2
% alpha is the degree of inversion imparted by the pulse

% T2 effects:

T2loss_bir8     = exp(-eTE*r2tis);
T2loss_art_bir8 = exp(-eTE*r2art) ;

T2loss_ftvsi     = exp(-eTE_ftvsi*r2tis);
T2loss_art_ftvsi = exp(-eTE_ftvsi*r2art);

T2loss_art_ftvss = exp(-eTE_ftvsi*r2art);

% Degree of inversion of each pulse
switch(label_type)
    case 'BIR8inv'
        % if the first pulse is  BIR8-inv : inversion efficiency with  T2 effects
        % tissue gets inverted too:  M(n) =  M(n-1) * (-alpha1_tis_sel);
        %
        alpha1_tis_sel  = -T2loss_bir8;
        alpha1_tis_ns   = -T2loss_bir8;
        alpha1_art_sel  = 0;
        alpha1_art_ns   = -T2loss_art_bir8;
        pulse1_dur    = 0.0272;
        %}
        
    case 'BIR8'
        % first pulse BIR8-sat : T2 effects
        % tissue does this: M(n) =  M(n-1) * alpha1_tis_sel
        %
        alpha1_tis_sel = T2loss_bir8;
        alpha1_tis_ns  = T2loss_bir8;
        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
        alpha1_art_sel = 0;
        alpha1_art_ns  = T2loss_art_bir8;
        pulse1_dur    = 0.0272;


        %
    case 'FTVSI-sinc'
        % first pulse BIR8-sat : T2 effects
        % tissue does this: M(n) =  M(n-1) * alpha1_tis_sel
        %
        alpha1_tis_sel = -T2loss_ftvsi;
        alpha1_tis_ns  = -T2loss_ftvsi;
        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
        alpha1_art_sel = T2loss_art_ftvsi;
        alpha1_art_ns  = -T2loss_art_ftvsi;

        pulse1_dur    = 0.0690 ;

    case 'FTVSS'
        % tissue gets saturated and arteries do not in the selective case
        % all gets saturates in the non-selective case
        alpha1_tis_sel = 0;
        alpha1_tis_ns  = 0;

        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
        alpha1_art_sel = T2loss_art_ftvss;
        alpha1_art_ns  = 0;
        pulse1_dur    = 17846 * 4e-6 ;

end
pulse2_dur = 0.0272;
pulse1_length = round(pulse1_dur/dt);
pulse2_length = round(pulse2_dur/dt);

% second pulse : T2 effects only on the tissue
alpha2_tis_sel = T2loss_bir8;
alpha2_tis_ns  = T2loss_bir8;

% second pulse : SATURATion or T2 effects only on the artery
alpha2_art_sel = 0;
alpha2_art_ns  = (T2loss_art_bir8);


% include B1 errors (bad efficiency) LHG 4/11/23
% alpha1_art_sel = alpha1_art_sel * b1err;
% alpha1_art_ns  = alpha1_art_ns * b1err;
% 
% alpha2_art_sel = alpha2_art_sel * b1err;
% alpha2_art_ns  = alpha2_art_ns * b1err;

% alpha1_tis_sel = alpha1_tis_sel * b1err;
% alpha1_tis_ns  = alpha1_tis_ns * b1err;
% 
% alpha2_tis_sel = alpha2_tis_sel * b1err;
% alpha2_tis_ns  = alpha2_tis_ns * b1err;

% bolus arrival time
nbat = round(bat/dt);

Nframes = length(del2);
obs = zeros(Nframes,1);

% the VS pulses are executed by the pulse sequence in a separate core.  
% Here, we include them at the end of the first delay (t_adjust)
% their effect is 'felt' at the end
% del1 = del1 + pulse1_dur;

TR = del1 + del2 + del3 + RO_time;

begTR = cumsum(TR);
begTR = [0; begTR];  % beginning of each TR period
begTR = begTR(1:end-1);
duration = sum(TR) ;
Npts = round((duration) / dt);
timevec = linspace(0, duration, Npts);

% indicator function for  labeling pulses
% whether the pulse is velocity selectivity is specified by "labelcontrol"
vsfun = zeros(Npts,1);

inds = round((begTR + del1) /dt);
inds = sort(inds);
inds = inds + pulse1_length;

% Notes
% in the scanner file:
% labelcontrol = 1 means vs pulse
% labelcontrol = 0  means non-vs pulse
% labelcontrol = -1  means No pulse at all

% in THIS code :
% vsfun(n) = 1 means vs pulse
% vsfun(n) = 0  means No pulse at all
% vsfun(n) = -1  means non-vs pulse

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
% in the scanner :
% doArtSuppression = 1 means AS pulse
% doArtSuppression = -1  means No pulse at all
% doArtSuppression = 0  means non AS pulse (T2 only)

asfun = zeros(Npts,1);

% indicator function for tissue T2 effects from arterial suppression:
% When we apply only the RF of arterial saturation pulses, we don't crush
% the arterial signal but we still get some T2 suppression on the tissue.
tsfun = zeros(Npts,1);


% the indices of when the AS pulses can be applied:
inds = round((begTR + del1  + del2) /dt);
% the AS pulse is inside the ArtSup delay core in the scanner
% but their effect is is felt at the end of the pulse
inds = inds + pulse2_length;

% velocity selective AS pulse: destroys arterial magnetization
tmp = zeros(size(doArtSup));
tmp(doArtSup ==1 ) = 1;  % both zeros and ones produce T2 effects
as_inds = inds .* tmp;
as_inds = inds(as_inds>0);
asfun(as_inds) = 1;

% non- selective AS pulse: T2 effects on arterial magnetization
tmp = zeros(size(doArtSup));
tmp(doArtSup ==0 ) = 1;  % both zeros and ones produce T2 effects
as_inds = inds .* tmp;
as_inds = inds(as_inds>0);
asfun(as_inds) = -1;

% make an indicator for the beginning of readout
aqfun = zeros(Npts,1);
inds = round((begTR + del1 +  del2 + del3) /dt);
aqfun(inds) = 1;

%------------
% LHG :2.5.22 - adapt the effect of readout for both FSE and GRE
% approximate the  effects of the read out on the Mz of tissue
% what is the "steady state" of the magnetization after the readout pulses.
if RO_type == 'FSE'
    flip = deg2rad(90);
end

flip = flip*b1err;  % this doesn't seem to matter

cosflip = cos(flip);  % pre-compute this once to save time
%-------------

img_time = round(RO_time / dt);  % how long it takes to collect image
slice_time = round(img_time/Nkz);


MM=2;

% Create arterial input functions:
%
% Each prep pulse affects a bouls of blood.
% Create the "Bolus" as rectangular input multiplied by a T1 relaxation
% function.
%
% The readout also affects the magnetization - it's another bolus
%
% Note: at the end of each bolus, we restore the magentization with fresh spins
% ie -  skip calculating the decay in the loop for that
% time point and reset the label to zero instead

bolus_duration = 3 ; % 2.5; % assumed bolus duration created by profile of labeling pulse
% now in units of discrete samples:
bolus_length = round(bolus_duration/dt);

tt = linspace(0,bolus_duration, bolus_length);
bolus = exp(-tt*r1a)';

tt_aq = linspace(0, RO_time, img_time);  % bolus duration from imaging
% in FSE readout,  assume that the first readout pulse saturates the arterial spins
% and the refocusing train doesn't affect the label significantly
AQbolus = exp(-tt_aq*r1a)';  
% in GRE readout, assume the readout has no effect on the label
if RO_type=='GRE'
    AQbolus(:) = 0;
end
% initialize the arterial concentration of "label" to zero
Ma = ones(Npts, 1);

%Ma(1:bolus_length)= 0.5 * bolus;  % everything starts with a BGS pulse

for n = MM:length(timevec)
   % 
    % Labeling pulses (vsfun) effect on the arterial contents:
    if vsfun(n-1) == 1  % gradients ON - velocity selective
        Ma(n )= Ma(n-1)*alpha1_art_sel;
    end    
    if vsfun(n-1) == -1 % gradients OFF - NOT velocity selective
        Ma(n) =  Ma(n-1)*alpha1_art_ns;
    end
   
    % ArtSup pulses cause saturation of the blood.
    % saturation recovery curve
    if asfun(n-1) == 1
        Ma(n)= alpha2_art_sel;
    end
    
    % Assume that Artsup control pulses don't do anything to  arterial
    % blood.  T2 effects on tissue spins only.
    %{
    % ArtSup control pulses cause T2 attenuation
    if asfun(n-1) == -1
        if Ma(n-1) > 0.5
            Ma(n+1 : n+bolus_length)= Ma(n-1)*alpha2_art_ns * bolus;
        else
            Ma(n+1 : n+bolus_length)= Ma(n-1)*(0.5+alpha2_art_ns) * bolus;
        end
    end
    %}

    
%----- test 8.2.23 (no effect from global saturation pulse on the arteries)
%
%----- end test code 8.2.23 
    % LHG: 2/6/2022  adapt for FSE vs. GRE 
    % Effect of the readout is that it perturbs magnetization in arteries
    % as well as the tissue.
    % Assume that the tipdown (90) saturates the arteries
    % but that the refocusers don't affect the Mz component.
    % After the readout, there is a global saturation pulse.
    % This resets the magnetization in all compartments to zero.
    if aqfun(n)==1
%         if n>img_time+1
%             if RO_type == 'FSE'
%                 % tipdown pulse saturates arteries
%                 Ma(n+1  : n+img_time)= 0.5 * AQbolus ;
%             end
%         end

%         % BGS pulse after readout also saturates the arteries
%         Ma(n+img_time+BGS0_length+1 : ...
%             n+img_time+BGS0_length + bolus_length)= 0.5 * bolus;
    end
    %}

    % Now implement Bloch equation
    dMa = (1-Ma(n))*r1a;
    Ma(n+1) = Ma(n) + dMa*dt;

end



%Ma = Ma2;
doDispersion = 1;
if doDispersion

    % arterial exchange compartment:  dispersion and delay
    %D = exp(-(tt-bat).^2/bat*100);

    % change to gamma function (1/31/22)
    D = tt.* exp(-(tt/bat*2).^2);
    D = D / sum(D);
    Ma_ex = conv(Ma, D);
    Ma_ex = Ma_ex(1:length(Ma));
else
    Ma_ex = Ma;
end
% scale to include T1 decay:
Ma_ex = Ma_ex/max(Ma_ex) ;%* max(Ma)*exp(-bat*r1a);

% convert from to concentration to Magnetization
%Ma_ex = 1-2*Ma_ex;
%Ma = 1-2*Ma;
%Ma = Ma_ex;


% Tissue Magnetization time courses
% Do NOT assume a BGS pulse at the very begining
M = ones(Npts, 1);

aq = 1;
Ma0 = Ma(1);
start_readout = 0;
kz = 0;
batcount = 1;
start_readout = inf;

%precalculate loop variables
sinflip = sin(flip);

% padding added on 2020-09-25 for patching up 'array index out of bounds'
% error
aqfun = [aqfun; zeros(round(del3(end)/dt),1)];
del3 = [del3;del3(end)];
aqfun_plot = aqfun;

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
    % dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n-1) - f * M(n-1)/lambda;
    dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n-1) - f * M(n-1)/lambda;
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
    if  aqfun(n)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments
        
        % LHG 4/17/23 - added the missing abs()       
        obs(aq) = abs( (1-cbva)*M(n-1)*sinflip + cbva*Ma(n-1) * sinflip); 
        
        obs_t(aq) = M(n-2) * sinflip;
        obs_a(aq) = Ma(n-2)* sinflip;
        
        aq = aq +1;
        start_readout = n;  % mark the begining of the acquisition
        
        % effect of the readout pulse
        M(n) = M(n-1) * cosflip;
        kz = 1;
        aqfun_plot(n  : n+img_time) = 1;

    end


    % LHG:2/6/22 - adapted for FSE and GRE cases 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % readout causes partial saturation of spins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if  (n == start_readout + slice_time*kz)
       M(n) = M(n) * cosflip;
       kz = kz + 1;
       if kz == Nkz +1
           kz = 1;
       end
     end

     % effect of global sat pulse (BGS) at end of readout: reset the magnetization
     if  (n == (start_readout + img_time + BGS0_length))
         M(n) = 0;
         %M(n) = b1err-1;
         start_readout = 0;
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


%%
if dofigs
    
    %some times they don't match up.  rounding error somewhere?
    ed = length(M);
    labelfun = vsfun;
    labelfun = labelfun(1:length(timevec));
    aqfun_plot = aqfun_plot(1:length(timevec));
    
    %figure()
    subplot(311)
    area(timevec(1:ed),aqfun_plot(1:ed),'FaceColor', 0.8*[1 1 1]);
    hold on
    stem(timevec(1:ed),labelfun(1:ed),'r')
    stem(timevec(1:ed),asfun(1:ed),'m');
    
    legend('Readout', 'Label', 'ArtSup')
    hold off
    axis tight
    
    subplot(312)
    area(timevec(1:ed),aqfun_plot(1:ed),'FaceColor', 0.8*[1 1 1])
    hold on
    plot(timevec(1:ed),M(1:ed),'b')
    plot(timevec(1:ed),Ma(1:ed),'g')
    plot(timevec(1:ed),Ma_ex(1:ed),'r')
    legend('ReadOut','Tis', 'Art', 'ArtX')
    axis([min(timevec) max(timevec) -1 1])
    
    hold off
  %  axis tight
    grid on
    
    subplot(313)
    plot(obs_a)
    hold on
    plot(obs_t)
    plot(obs)
    %plot(cbva*obs_a(2:end))
    legend('Arterial', 'Tissue', 'Both') % , 'Art')
    hold off
    grid on
    
    drawnow()
end



return
