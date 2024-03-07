function  [obs Ma] = gen_signals_vs_230918(parms, aq_parms, dofigs, doSub, dt)
% this version is modified to agree with the timings of the new pulse
% sequence asl3dflex
% it also considers the effect of multiple pulses on the arterial
% compartment - splits into multiple populations
% this function returns the arterial Magnetization function Ma (before
% dispersion)

% constants
if nargin<5
    dt = 1e-3;      % seconds
    %dt = 5e-3;
    %dt = 0.01;      % seconds
end

r1a = 1/1.67;
lambda = 0.9;
tissDensity=1.05;
r2art  = 1/0.150 ;

eTE_ftvsi = 0.0313;     % effective echo time of FT-VSI sinc-pulse
eTE = 0.0207;           % effective echo time of sBIR8-VSS pulse

BGS0_length = round(2.5e-3 /dt); % (s)  duration of the global saturation pulse - magnetization reset.
BGS0_length = round(6.4e-3 /dt); % (s)  duration of the global saturation pulse - magnetization reset.
BGS0_length = 0;   % did not use bulk saturation

if (nargin ==0)
    % for testing purposes, here are some default tissue parameters:
    % Tissue parms:
    f=          0.01 ;
    cbva =      0.02;
    bat =       0.1 ;
    r1tis =     1/1.4  ;
    Mtis0 =     1 ;
    flip =      50*pi/180 ; % flip angle in radians
    r2tis =     1/0.090 ;

    % include a B1 error term in the labeling pulses
    % 0 error means multiplying by 1.
    % 1% error mean multiplying by 1.01
    b1err =     0.05;

    % processing defaults
    doSub = 1;
    dofigs = 1;

    % pulse sequence parameters
    Nframes = 8;
    label_type =    'BIR8inv'; %'FTVSS'; %'FTVSI-sinc'; % 'BIR8inv'; % 'BIR8'
    RO_type =       'GRE'; %'FSE';   % 'GRE'
    t_tags =        0;% 0.1*ones(Nframes,1);
    del1 =          2*ones(Nframes, 1);
    del2 =          1.2*ones(Nframes,1);
    del3 =          0.15*ones(Nframes,1);  % delay between AS pulse and acqusition
    labelcontrol =  zeros(Nframes,1);
    labelcontrol(2:2:end)= 1;
    labelcontrol(1:2) =   -1;
    order =         1;
    doArtSup =      ones(Nframes,1);
    %doArtSup(end/2+1:end)= 0;
    %doArtSup(1:2)=    0;
    % test case
    %{
    del1 = 1.5*ones(Nframes, 1);
    del2 = 0.5*ones(Nframes,1);
    %}
    Ma = [];

    Nkz = 1;
   RO_time = 0.03*Nkz;
 
else
    Nkz = 1;

    % pulse sequence parms:
%    t_tags =        aq_parms.t_tags;  % zero if you have a single pulse
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
    flip =          aq_parms.flip;

    if isfield(aq_parms, 'Ma')  % maybe the input function is already precomputed.
        Ma = aq_parms.Ma;
    else
        Ma=[];
    end

    % Tissue parms
    f = parms.f;
    Mtis0 = parms.Mtis0;
    cbva = parms. cbva;
    bat =  parms.bat;
    r1tis =  parms.r1tis;
    r2tis = parms.r2tis;

    % include a B1 error term in the labeling pulses
    % 0 error means multiplying by 1.
    % 1% error mean multiplying by 1.01
    b1err = parms.b1err;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Translate labeling error (B1err) into a scaling factor
b1err= b1err +1;

% Apply flip angle error factor
flip = flip*b1err;  

% check the order.  The 'order' flag is a quick way to reverse the order
% ie - VS becomes non-VS and non-VS becomes VS
tmplabelcontol= labelcontrol;
if order==2  %  2023.08.08....original should be  2 (just checking to see if there is an error)
    tmplabelcontrol(labelcontrol==0) = 1;
    tmplabelcontrol(labelcontrol==1) = 0;
    tmplabelcontrol(labelcontrol==-1) = -1;

    labelcontrol = tmplabelcontrol;
end


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
        % tissue gets inverted too:  M(n) =  M(n-1) * alpha1_tis_sel;
        %
        alpha1_tis_sel  = -T2loss_bir8;
        alpha1_tis_ns   = -T2loss_bir8;

        % artery input function does this:  Ma(n) =  alpha1_art_sel * bolus;
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

pulse2_dur = 0.0272;  % always the same - BIR8

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


% adjust the delay times from the file to account for the time it takes to
% play the prep pulse, then we play the delay.
% we'll make the effects of the prep pulse happen at pulse1_dur
de1 = del1 + BGS0_length*dt;
del2 = del2 + pulse1_dur;
del3 = del3 + pulse2_dur;

% No BGS) pulse


TR = del1 + del2 + del3 + RO_time;

begTR = cumsum(TR);
begTR = [0; begTR];  % beginning of each TR period
begTR = begTR(1:end-1);
duration = sum(TR) ;
Npts = round((duration) / dt);
timevec = linspace(0, duration, Npts);

% indicator function for  labeling pulses
% whether the pulse is velocity selective is specified by "labelcontrol"
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

vs_inds = inds(:) .* tmp(:);  % vel. selective pulse
vs_inds = inds(vs_inds>0);

tmp = zeros(size(labelcontrol));
tmp(labelcontrol==0) = 1;
ns_inds = inds(:) .* tmp(:);  % non-selective pulse
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

% non- selective AS pulse: T2 effects on magnetization
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

cosflip = cos(flip);  % pre-compute this once to save time
%-------------

img_time = round(RO_time / dt);  % how long it takes to collect image
slice_time = round(img_time/Nkz);



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
bolus_duration = 2.5; % assumed bolus duration created by profile of labeling pulse
% now in units of discrete samples:
bolus_length = round(bolus_duration/dt);
tt = linspace(0,bolus_duration, bolus_length);
bolus = exp(-tt*r1a)';


% in FSE readout,  assume that the first readout pulse saturates the arterial spins
% but the refocusing train doesn't affect the label significantly
tt_aq = linspace(0, RO_time, img_time);  % bolus duration produced by imaging tip pulses
AQbolus = exp(-tt_aq*r1a)';
% in GRE readout, assume the readout has no effect on the label
if RO_type=='GRE'
    AQbolus(:) = 0;
end


% Note that the arterial input function could be precomputed already
% if it's empty, then we'll need to compute it now;
if isempty(Ma)
    % initialize the arterial concentration of "label" to one
    Ma = ones(Npts, 1);
    Ma_fresh = ones(Npts, 1);  % these are the fresh arterial spins
    % that haven't seen any pulses yes

    %Ma(1:bolus_length)= 0.5 * bolus;  % everything starts with a BGS pulse

    last_event = 1;
    flow_time = 0;
    kz=0;
    start_readout=inf;

    for n = 2:length(timevec)
        %
        % Labeling pulses (vsfun) effect on the arterial contents:
        if vsfun(n-1) == 1  % gradients ON - velocity selective
            Ma(n )= Ma(n-1)*alpha1_art_sel;
            % we now combine the two populations according to how much time
            % elapsed since the last pulse.
            % the weight is a sigmoid centered around the bolus duration
            %Ma_fresh(n) = alpha1_art_sel;
            %flow_time = dt*(n-last_event);
            %last_event = n;
            %W = 1/(1 + exp(-(3.0-flow_time)/0.1));
            %Ma(n) = W*Ma(n) + (1-W)*Ma_fresh(n);

        end
        if vsfun(n-1) == -1 % gradients OFF - NOT velocity selective
            Ma(n) =  Ma(n-1)*alpha1_art_ns;
            % we now combine the two populations according to how much time
            % elapsed since the last pulse.
            % the weight is a sigmoid centered around the bolus duration
%             Ma_fresh(n) =  alpha1_art_ns;
%             flow_time = dt*(n-last_event);
%             last_event = n;
%             W = 1/(1 + exp(-(3.0-flow_time)/0.1));
%             Ma(n) = W*Ma(n) + (1-W)*Ma_fresh(n);

        end
        % ArtSup pulses cause saturation of the blood.
        % saturation recovery curve
        if asfun(n-1) == 1
            Ma(n)= alpha2_art_sel;
            % A thought?: 
            % we now combine the two populations according to how much time
            % elapsed since the last pulse.
            % the weight is a sigmoid centered around the bolus duration
            %Ma_fresh(n)= alpha2_art_sel;
            %flow_time = dt*(n-last_event);
            %last_event = n;
            %W = 1/(1 + exp(-(3.0-flow_time)/0.1));
            %Ma(n) = W*Ma(n) + (1-W)*Ma_fresh(n);

        end

        % Assume that Artsup control pulses don't do anything to  arterial
        % blood.  T2 effects on tissue spins only.
        %
    % ArtSup control pulses cause T2 attenuation
    if asfun(n-1) == -1
            Ma(n) =  Ma(n-1)*alpha1_art_ns;
    end
        %}


        %----- test 8.2.23 (no effect from global saturation pulse on the arteries)
        %
        %----- end test code 8.2.23
        % LHG: 2/20/2024  adapt for FSE vs. GRE
        % Effect of the readout is that it perturbs magnetization in arteries
        % as well as the tissue.
        if aqfun(n)==1
            start_readout = n;  % mark the begining of the acquisition

            % effect of the readout pulse
            Ma(n) = Ma(n-1) * cosflip;
        end
        if  (n == start_readout + slice_time*kz)
            Ma(n) = Ma(n) * cosflip;
            kz = kz + 1;
        end
        if kz==Nkz
            kz=0;
        end

        % BGS pulse after readout also saturates the arteries
        %Ma(n+img_time+BGS0_length+1 : ...
        %    n+img_time+BGS0_length + bolus_length)= 0.5 * bolus;


        % Now implement Bloch equation
        dMa = (1-Ma(n))*r1a;
        Ma(n+1) = Ma(n) + dMa*dt;

        %dMa_fresh = (1-Ma_fresh(n))*r1a;
        %Ma_fresh(n+1) = Ma_fresh(n) + dMa_fresh*dt;


    end
end

% scale to match the spin density
Ma = Ma*Mtis0;

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
%Ma_ex = Ma_ex/max(Ma_ex) ;%* max(Ma)*exp(-bat*r1a);



% Tissue Magnetization time courses
% Do NOT assume a BGS pulse at the very begining
M = ones(Npts, 1);

aq = 1;
Ma0 = Ma(1);
start_readout = inf;
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

MM=2; %0.5/dt;

for n = nbat+MM+1:length(timevec)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The diff eqs governing the three pools of protons
    % in the voxel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %  blood in artery
    %dMa = (Ma0 - Ma(n-1))*r1a ;
    %Ma(n) = Ma(n-1) + dMa*dt;

    % the modified Bloch equation has
    % t1 decay, , inflow, outflow
   
    % dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n-1) - f * M(n-1)/lambda;
    dM = (Mtis0 - M(n-1))*r1tis  + f * Ma_ex(n-1);

    M(n) = M(n-1) + dM*dt;

    %%%%%%%%%%%%%%%%%%%
    % update the Tissue and Arterial Compartments when  pulses are
    % applied
    %%%%%%%%%%%%%%%%%%%

    % selective pulse
    if vsfun(n-1) == 1
        M(n-pulse1_length:n) =  M(n-pulse1_length) * alpha1_tis_sel;
    end
    % non-selective pulse
    if vsfun(n-1) == -1
        M(n-pulse1_length:n) =  M(n-pulse1_length) * alpha1_tis_ns;
    end

    % T2 effects from the arterial selective saturation pulse
    % whether the crusher gradients are on or not.
    if asfun(n-1) ~= 0
        M(n - pulse2_length:n) =  M(n-pulse2_length) * (alpha2_tis_ns); % T2 effect only : no inversion
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % Readout  from the different pools
    % dependence on flip angle
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if  aqfun(n)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments

        % LHG 4/17/23 - added the missing abs()
        obs(aq) = abs( (1-cbva) * tissDensity * M(n-1)*sinflip + cbva*Ma(n-1) * sinflip);

        obs_t(aq) = M(n-1) * sinflip;
        obs_a(aq) = Ma(n-1)* sinflip;
        
        aq = aq +1;
        start_readout = n;  % mark the begining of the acquisition

        % effect of the readout pulse
        M(n) = M(n-1) * cosflip;
        %kz = 1;
        aqfun_plot(n  : n+img_time) = 1;

    end


    % LHG:2/6/22 - adapted for FSE and GRE cases
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % readout causes partial saturation of spins
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  (n == start_readout + slice_time*kz)
        M(n) = M(n) * cosflip;
        kz = kz + 1;
    end
    if kz == Nkz
        kz = 0;
    end

    %{
    % effect of global sat pulse (BGS) at end of readout: reset the magnetization
    if  (n == (start_readout + img_time + BGS0_length))
        M(n) = 0;
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
