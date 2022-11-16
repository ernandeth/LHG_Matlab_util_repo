function  obs = gen_signals_ps_v2(parms, tags, timing_parms, dofigs,doSub)

% constants
global dt
dt = 0.005;      % seconds

% dt = 0.01;      % seconds
% r1a = 1/1.67;
r1a = 1/1.7;
%r1a = 0.00001;

lambda = 0.9;
aqwindow = 0.034 ;
alpha = 0.85;

if nargin==0
    % for testing purposes, here are some default parameters:
    f=         60 /6000;
    mtis0 =     1;
    cbva =      0.01 ;
    eta=   0.8 ;
    bat_tis = 1.2;
    bat_art = 0.6;
    r1tis =     1/1.4  ;
    flip =      70*pi/180 ; % flip angle in radians
    Disp =      40;
    r1blood = 1/1.7;
    
    Nframes = 8;
    aqwindow = aqwindow * ones(Nframes, 1);
    % timing for straight up  PCASL :
    t_delay  = 0.2 * ones(Nframes,1);
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
    eta= parms.eta;
    bat_tis = parms.bat_tis;
    bat_art = parms.bat_art;
    r1tis = parms.r1tis;
    flip = parms.flip;
    Disp = parms.Disp;
    r1blood = parms.r1blood;
    
    t_delay =           timing_parms.t_delay;
%     t_tag =              timing_parms.t_tag;
    t_tag = tags;
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

% changed to zeros 11 July
M_tis = ones(Npts, 1); % why ones? ask!
M_blood = M_tis;

labelfun = zeros(Npts,1);
aqfun = zeros(Npts,1);
obs = zeros(Nframes,1);
%%
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

% arterial input magnetization
inMa_tis = labelfun;
inMa_tis(inMa_tis<0) = 0;



% Arterial dispersion and decay Kernel from labeling plane to voxel
% arterial compartment only
ttmp = linspace(0,10, 10/dt)';
% decay = ones(size(ttmp));
decay = exp(-ttmp*r1a);

doDispersion=0;

if doDispersion
    art_kernel = ttmp.^(bat*Disp) .* exp(-Disp *ttmp);
    
    art_kernel(ttmp <= 0)=0;
    art_kernel(isnan(art_kernel)) = 0;
    art_kernel = abs(art_kernel);
    art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
    art_kernel = art_kernel .* decay;
    art_kernel(isnan(art_kernel)) = 0;
    
    % second delay: time spent in the arteries before exchange:
    % If no dispersion in kernel, only T1 decay and a time shift
    art_kernel2 = ttmp.^(bat_tis*Disp) .* exp(-Disp *ttmp);
    art_kernel2(isnan(art_kernel2)) = 0;
    art_kernel2(ttmp <= 0)=0;
    art_kernel2 = abs(art_kernel2);
    art_kernel2 = art_kernel2 / sum(art_kernel2);  %normalize the mass in the dispersion
    art_kernel2 = art_kernel2 .* decay;
    art_kernel2(isnan(art_kernel2)) = 0;
else
    
    %------------  LHG: 1.21.16 ------
    % If no dispersion in kernel, only T1 decay and a time shift
    %
    % delay to the local arteries from labeling plane
    % (note the +1 shift allows for eta, bat2 ==0 ) 
    
    % for the output of pass through artery
    art_kernel(:) =  zeros(size(ttmp));
    art_kernel(round((bat_art+eta)/dt) + 1) = exp(-(bat_art+eta)*r1a);
    
    % for the input to the tissue compartment
    art_kernel2 = art_kernel;
    art_kernel2(:) = 0;
    art_kernel2(round( bat_tis /dt ) +1) = exp(-bat_tis*r1a);
    
    % for the input to the pass through artery
    art_kernel3 = art_kernel;
    art_kernel3(:) = 0;
    art_kernel3(round( bat_art /dt ) +1) = exp(-bat_art*r1a);
end

% -----------------

% arterial CBV magnetization after bolus arrival and dispersion.
M_out = labelfun;
M_out(M_out<0) = 0;
inMa_tis = M_out;
inMa_art = M_out;
% label content in the arterial compartment
M_out =  conv(M_out, art_kernel);
M_out = M_out(1:Npts);

% label input to the tissue (after the bolus arrival time)
inMa_tis = conv(inMa_tis, art_kernel2);
inMa_tis = inMa_tis(1:Npts);

inMa_art = conv(inMa_art,art_kernel3);
inMa_art = inMa_art(1:Npts);

% Z magnetization in the arteries 
M_out =  mtis0*(1 - 2*alpha*M_out);         % output of the arteries
inMa_tis = mtis0*(1 - 2*alpha*inMa_tis);     % goes into the tissue
inMa_art = mtis0*(1 - 2*alpha*inMa_art);    % input to the arteries

aq = 1;

%%
for n=2:length(timevec)
        
    % the modified Bloch equation has
    % t1 decay,  magnetization transfer, inflow, outflow

    % Equation [2] Su et al
    dM_tis = (mtis0 - M_tis(n-1))*r1tis + f * inMa_tis(n-1)  - f * M_tis(n-1)/lambda;
    M_tis(n) = M_tis(n-1) + dM_tis*dt;
    
    % Equation [3] Su et al
    dM_blood = (mtis0 - M_blood(n-1))*r1blood + (inMa_art(n-1)-M_out(n-1))/eta;
    M_blood(n) = M_blood(n-1) + dM_blood*dt;
    
    
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if aqfun(n-1)==1
        % the observed signal on the xy plane:  Tissue + Blood compartments 
        obs(aq) = ( ( (1-cbva)*M_tis(n-1)*1.05 ) + ( cbva*M_blood(n-1) ) ) * sin(flip);
         
        aq = aq+1;
        
        % what is left on the z axis after the rf pulse.
        M_tis(n) =  M_tis(n-1)  * cos(flip);
        M_out(n) =  M_out(n-1)  * cos(flip);
        inMa_tis(n) =  inMa_tis(n-1)  * cos(flip);
        
    end
    
end

if dofigs
    

    labelfun = labelfun(1:length(timevec));
    aqfun = aqfun(1:length(timevec));
    
    figure(10)
    subplot(311)
    area(timevec,labelfun)
    hold on
    stem(timevec,aqfun,'r');
    hold off
    
    subplot(312)
    plot(timevec,M_tis,'b')
    hold on
    plot(timevec,M_blood,'r')
    hold on
    
    plot(timevec,M_out,'--b')
    plot(timevec,inMa_tis,'--g')
    plot(timevec,inMa_art,'--r')
    grid on
    hold off
    axis tight
    
    subplot(313)
    plot(obs)
    
end

if doSub
    obs = obs(2:2:end) - obs(1:2:end);
end


return
