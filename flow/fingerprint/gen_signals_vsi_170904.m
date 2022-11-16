function  obs = gen_signals_vsi_170904(parms, timing_parms, dofigs,doSub)

% constants
dt = 0.001;      % seconds
dt = 0.0001;      % seconds

% dt = 0.01;      % seconds
% r1a = 1/1.67;
r1a = 1/1.7;
%r1a = 0.00001;

lambda = 0.9;
aqwindow = 0.034 ;
alpha = 0.8;
vsi = 0;

if nargin==0
    % for testing purposes, here are some default parameters:
    f=         60 /6000;
    mtis0 =     1 ;
    cbva =      0.01 ;
    bat =       0.3 ;
    bat2 =      1.2;
    kfor =      0 ;
    r1tis =     1/1.4  ;
    flip =      60*pi/180 ; % flip angle in radians
    Disp =      40;
    TR = 4;
    
    Nframes = 10;
    aqwindow = aqwindow * ones(Nframes, 1);
    % timing for straight up  VSASL :
    t_delay  = 0.35 * ones(Nframes,1);            % time between arterial crusher and acquisition
    t_tag = 1.5*ones(Nframes,1);                  % time between label pulse and arterial crusher
    t_adjust = 1.5*ones(Nframes,1) - aqwindow;    % time between acquisition and label pulse
    t_tag = 1.5*randn(Nframes,1);                 % time between label pulse and arterial crusher
    
    order = 1;  % control-tag
    Nlabel_group = 1;  % Number of tags in a row  (added this on 5/7/15)
    
    tmp = ones(size(t_adjust));
    tmp(1:2:end) = 0;
    isLabel = tmp';
    
    doSub = 1;
    dofigs = 1;
    
else
    
    % get parms from the input structure
    mtis0 = parms.mtis0;
    f = parms.f;
    cbva = parms.cbva;
    bat = parms.bat;
    bat2 = parms.bat2;
    kfor = parms.kfor;
    r1tis = parms.r1tis;
    flip = parms.flip;
    Disp = parms.Disp;
    
    
    t_delay =           timing_parms.t_delay;
    t_tag =              timing_parms.t_tag;
    t_adjust =          timing_parms.t_adjust;
    Nlabel_group =  timing_parms.Nlabel_group;
    order =               timing_parms.order;
    isLabel =            timing_parms.isLabel ;
    aqwindow =       timing_parms.t_aq;
end

if order==2
    isLabel= ~isLabel;
end

Nframes = length(t_delay);

TR = t_adjust + t_tag + t_delay + aqwindow;

duration = sum(TR) +1;
% the extra  second is to prevent the last sample from getting chopped off.

Npts = floor(duration / dt);
timevec = linspace(0, duration, Npts);

M = ones(Npts, 1);
Ma = ones(Npts, 1);
inMa = ones(Npts, 1);

vsifun = zeros(Npts,1);
aqfun = zeros(Npts,1);
vssfun = zeros(Npts,1);
obs = zeros(Nframes,1);
satfun = zeros(Npts,1); % this captures the saturation 
                        % from sampling the signal during acquisition

                        
begTR=[0 ; cumsum(TR(1:end-1))];
% times for VSI pulse
inds =  floor((begTR + t_adjust) /dt);
vsifun(inds)=1;
% times for NonVSI pulse
vsifun(inds(isLabel==0)) = -1;

% times for the end of the input function (VSS pulse)
inds = floor((begTR + t_adjust + t_tag) /dt);
vssfun(inds) = 1;

% times for acquisition
inds = floor((begTR + t_adjust + t_tag + t_delay) /dt);
aqfun(inds) = 1;
    
% times for saturation of stationary spins by the acquisition window
aqpts = floor(aqwindow/dt);
for n=0:aqpts
    inds = floor((begTR + t_adjust + t_tag + t_delay ) /dt + n);
    satfun(inds) = 1;
end

aq=1;
Ma0 = 1;
for n=2:length(timevec)
    
    % arterial blood in voxel.
    dMa = (Ma0 - Ma(n-1))*r1a;
    Ma(n) = Ma(n-1) + dMa*dt;
    
    % arterial blood going into tissue after delay
    dMa = (Ma0 - inMa(n-1))*r1a;
    inMa(n) = inMa(n-1) + dMa*dt;
    
    % the modified Bloch equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dM = (mtis0 - M(n-1))*r1tis  - kfor*M(n-1) + f * inMa(n-1)  - f * M(n-1)/lambda;
    M(n) = M(n-1) + dM*dt;
    
    % velocity selective inversion pulse
    if vsifun(n-1) == 1
        M(n) = -alpha*M(n-1);
        Ma(n) = alpha*Ma(n-1);
    end
    
    % non-selective inversion pulse
    if vsifun(n-1) == -1
        Ma(n) = -alpha*Ma(n-1);
        M(n) = -alpha*M(n-1);
        vsi=1;
    end
    
    % sat pulse to determine the end of the input function
    if vssfun(n-1)==1
        Ma(n) = 0;
        vsi = 1;
    end
    
    % after BAT, the contents of the arterial compartment (Ma)
    % become the input funtion (inMa)
    nbat = round(bat/dt);
    if n > nbat
        % velocity selective
        if vsifun(n-nbat) == -1
            inMa(n) = Ma(n);
        end
        
        % velocity non-selective
        if vsifun(n-nbat) == 1
            inMa(n) = Ma(n);
        end
        
        % sat. pulse
        if vssfun(n-nbat)==1 
            inMa(n) = Ma(n);
        end
    end
    
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if  aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments
        obs(aq) = ((1-cbva)*M(n-1) + cbva*Ma(n-1)) * sin(flip);
        
        aq = aq+1;
    end
    
    if  satfun(n-1)==1

        % what is left on the z axis after sampling the signal with a 3D acquisition
        M(n) =  M(n-1)*cos(flip);
        Ma(n) =  Ma(n-1)*cos(flip);
        inMa(n) =  inMa(n-1)*cos(flip);

    end
    
end

if doSub
    obs = obs(2:2:end) - obs(1:2:end);
end

if dofigs
    
    %some times they don't match up.  rounding error somewhere?
    
    labelfun = vsifun;
    labelfun = labelfun(1:length(timevec));
    aqfun = aqfun(1:length(timevec));
    
    figure(1)
    subplot(311)
    stem(timevec,labelfun)
    hold on
    stem(timevec,aqfun,'r');
    stem(timevec,vssfun,'k');
    legend('Label','Aq','vasc. sup')
    hold off
        axis tight

    subplot(312)
    plot(timevec,M)
    hold on
    plot(timevec,Ma,'r')
    plot(timevec,inMa,'g')
    legend('Tis', 'Art', 'Art Input')
    hold off
    axis tight
    grid on
    
    subplot(313)
    plot(obs)
    grid on
end



return
