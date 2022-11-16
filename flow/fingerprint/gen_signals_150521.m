function  obs = gen_signals_150521(parms, timing_parms, dofigs,doSub)

% constants
dt = 0.001;      % seconds

dt = 0.005;      % seconds
r1a = 1/1.67;

lambda = 0.9;
aqwindow = 0.034 ;
alpha = 0.85;

if nargin==0
    % for testing purposes, here are some default parameters:
    f=         60 /6000;
    mtis0 =     1 ;
    cbva =      0.02 ;
    bat =   0.8 ;
    bat2 = 0.5;
    kfor =      1e-2 ;
    r1tis =     1/1.4  ;
    flip =      90*pi/180 ; % flip angle in radians
    Disp =      40;
    
    Nframes = 8;
    aqwindow = aqwindow * ones(Nframes, 1);
    % timing for straight up  PCASL :
    t_delay  = 1.6 * ones(Nframes,1);
    t_tag = 2*ones(Nframes,1);
    t_adjust = 0.4*ones(Nframes,1) - aqwindow;
    order = 1;  % control-tag
    Nlabel_group = 1;  % Number of tags in a row  (added this on 5/7/15)
    
    tmp = ones(size(t_adjust));
    tmp(1:2:end) = 0;
    isLabel = tmp';
       
    doSub = 0;
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

TR = t_adjust + t_delay  + t_tag + aqwindow;

duration = sum(TR) +1; 
% the extra  second is to prevent the last sample from getting chopped off.

Npts = floor(duration / dt);
timevec = linspace(0, duration, Npts);

M = ones(Npts, 1);
Ma = ones(Npts, 1);

labelfun = zeros(Npts,1);
aqfun = zeros(Npts,1);
MTfun = zeros(Npts,1);
obs = zeros(Nframes,1);

n1 = 1;
for n=1:Nframes
    n2 = n1 +  round( t_adjust(n)/dt)  ; % labeling begins
    n3 = n2 + round( t_tag(n)/dt )    ;  % labeling ends
    n4 = n3 + round(t_delay(n)/dt )   ; % flip angle happens
    n5 = n4 + round( aqwindow(n)/dt )   ; % end of the aquisition
    
    % [n n2 n3 t_delay(n) t_tag(n) t_adjust(n)]
    
    if isLabel(n)==1
        labelfun(n2:n3) = 1;
    end
    if isLabel(n)==0
        labelfun(n2:n3) = -1;
    end

    aqfun(n4) = 1;
    
    n1 = n1 + round(TR(n)/dt);  % begin the next frame
end


% MT happens in both control and labeling cases
kfor = kfor * abs(labelfun);

% arterial input magnetization
inMa = labelfun;
inMa(inMa<0) = 0;

% arterial CBV magnetization
Ma = labelfun;
Ma(Ma<0) = 0;

% Arterial dispersion and decay Kernel from labeling plane to voxel.
ttmp = linspace(0,10, 10/dt)';
decay = ones(size(ttmp));
decay = exp(-ttmp*r1a);

art_kernel = ttmp.^(bat*Disp) .* exp(-Disp *ttmp);
art_kernel(isnan(art_kernel)) = 0;
art_kernel(ttmp <= 0)=0;
art_kernel = abs(art_kernel);
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
art_kernel = art_kernel .* decay;

%
 %------------  LHG: 1.21.16 ------
% If no dispersion in kernel, only T1 decay and a time shift
%{
% delay to the local arteries from labeling plane
art_kernel(:) = 0 ;
art_kernel(round(bat/dt)) = exp(-bat*r1a);
%}

% second delay: time spent at the arteries before exchange:  local transit time
% this kernel is a 200 ms. delay and a little more T1 relaxation.
art_kernel2 = art_kernel;
art_kernel2(:) = 0;
art_kernel2(round( bat2 /dt )) = exp(-bat2*r1a);

%}
% -----------------

% label content in the arterial compartment
Ma =  conv(Ma, art_kernel);
Ma = Ma(1:Npts);

% label input to the tissue (after the local transit time)
inMa = conv(Ma, art_kernel2);
inMa = inMa(1:Npts);

% Z magnetization in the arteries 
Ma = mtis0- 2*alpha*Ma;
inMa = mtis0- 2*alpha*inMa;

aq = 1;
for n=2:length(timevec)
        
    % the modified Bloch equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dM = (mtis0 - M(n-1))*r1tis  - kfor(n-1)*M(n-1) + f * inMa(n-1)  - f*M(n-1)/lambda;
    M(n) = M(n-1) + dM*dt;
    
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments 
        obs(aq) = ((1-cbva)*M(n-1) + cbva*Ma(n-1)) * sin(flip);
         
        aq = aq+1;
        
        % what is left on the z axis after the rf pulse.
        M(n) =  M(n-1)  * cos(flip);
        Ma(n) =  Ma(n-1)  * cos(flip);
        
    end
    
end

if dofigs
    
    %some times they don't match up.  rounding error somewhere?
    labelfun = labelfun(1:length(timevec));
    aqfun = aqfun(1:length(timevec));

    figure(1)
    subplot(311)
    area(timevec,labelfun)
    hold on
    stem(timevec,aqfun,'r');
    hold off
    
    subplot(312)
    plot(timevec,M)
    hold on
    plot(timevec,inMa,'r')
    hold off
    axis tight
    
    subplot(313)
    plot(obs)
    
end

if doSub
    obs = obs(2:2:end) - obs(1:2:end);
end


return
