function  parms_result = fit_signals_160608( timing_parms, dofigs, data);

% constants
dt = 0.001;      % seconds
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
    parms.f=         60 /6000;
    parms.mtis0 =     1 ;
    parms.cbva =      0.01 ;
    parms.bat =   1.2 ;
    parms.bat2 = 1.8;
    parms.kfor =      1e-2 ;
    parms.r1tis =     1/1.4  ;
    parms.flip =      70*pi/180 ; % flip angle in radians
    parms.Disp =      200;
    
    Nframes = 50;
    aqwindow = aqwindow * ones(Nframes, 1);
    % timing for straight up  PCASL :
    t_delay  = 0.05 * ones(Nframes,1);
    t_tag = linspace(0.05, 2, Nframes);
    t_tag = t_tag(randperm(Nframes))';
    t_adjust = 0.05*ones(Nframes,1) - aqwindow;
    order = 1;  % control-tag
    Nlabel_group = 1;  % Number of tags in a row  (added this on 5/7/15)
    
%     tmp = ones(size(t_adjust));
%     tmp(1:2:end) = 0;
    tmp = round(rand(size(t_adjust)));
    isLabel = tmp';
    
    doSub = 0;
    dofigs = 1;
    
else
    t_delay =           timing_parms.t_delay;
    t_tag =             timing_parms.t_tag;
    t_adjust =          timing_parms.t_adjust;
    Nlabel_group =      timing_parms.Nlabel_group;
    order =             timing_parms.order;
    isLabel =           timing_parms.isLabel ;
    aqwindow =          timing_parms.t_aq;
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

experiment.labelfun = labelfun;
experiment.aqfun = aqfun;

parms_vector =  [ 
    parms.f ;
%    parms.mtis0 ;
    parms.cbva ;
    parms.bat  ;
    parms.bat2 ;
%    parms.kfor ;
    parms.r1tis ;
    parms.flip ; 
%    parms.Disp ;
    ];

LB = [0 ;  0; 0.2; 0.2; 0.35; deg2rad(20)];
UB = [0.02 ; 0.05; 3; 3; 2; deg2rad(70)];
initguess = parms_vector * 0.8;

hold off

% Make some fake data
obs = signal_cost(parms_vector, experiment,  []);
plot(obs)
obs = obs + 0.0000005*randn(size(obs));
hold on; plot(obs, '*'); 

% try to fit it now
myfun = @(parms_vector) signal_cost( parms_vector, experiment,  obs);

opts = optimset('lsqnonlin');

parms_result = lsqnonlin(myfun, initguess,  [], [], opts) ;

% alt. syntax?
% parms_result = lsqnonlin('signal_cost', initguess,  LB, UB, opts,...
%     experiment,  obs) ;

predicted = signal_cost(parms_result, experiment,  []);
plot(predicted, 'k');


parms_vector

return

%%
function result = signal_cost(parms_vector, experiment, data)
dofigs = 0;
% constants
dt = 0.001;      % seconds
dt = 0.005;      % seconds

% dt = 0.01;      % seconds
% r1a = 1/1.67;
r1a = 1/1.7;
%r1a = 0.00001;

labelfun = experiment.labelfun;
aqfun = experiment.aqfun ;

Npts = length(labelfun);
duration = Npts * dt;

lambda = 0.9;
aqwindow = 0.034 ;
alpha = 0.85;

% extract parms from the input structure
f = parms_vector(1);
cbva = parms_vector(2);
bat  = parms_vector(3);
bat2 = parms_vector(4);
r1tis = parms_vector(5);
flip = parms_vector(6);

Disp =40;
mtis0 = 1;
kfor =0;

% MT happens in both control and labeling cases
kfor = kfor * abs(labelfun);

% arterial input magnetization
inMa = labelfun;
inMa(inMa<0) = 0;

% Tissue Magnetization
M = ones(Npts, 1);

% Arterial compartment (flow-thru) magnetization
Ma = ones(Npts, 1);

% Arterial dispersion and decay Kernel from labeling plane to voxel
% arterial compartment only
ttmp = linspace(0,10, 10/dt)';
% decay = ones(size(ttmp));
decay = exp(-ttmp*r1a);


art_kernel = ttmp.^(bat*Disp) .* exp(-Disp *ttmp);

art_kernel(ttmp <= 0)=0;
art_kernel(isnan(art_kernel)) = 0;
art_kernel = abs(art_kernel);
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
art_kernel = art_kernel .* decay;
art_kernel(isnan(art_kernel)) = 0;

% Arterial dispersion and decay Kernel from labeling plane to voxel
% tissue/capillary compartment only
art_kernel2 = ttmp.^(bat2*Disp) .* exp(-Disp *ttmp);

art_kernel2(ttmp <= 0)=0;
art_kernel2 = abs(art_kernel2);
art_kernel2 = art_kernel2 / sum(art_kernel2);  %normalize the mass in the dispersion
art_kernel2 = art_kernel2 .* decay;
art_kernel2(isnan(art_kernel2)) = 0;


% arterial CBV magnetization after bolus arrival and dispersion.
Ma = labelfun;
Ma(Ma<0) = 0;
inMa = Ma;

% label content in the arterial compartment
Ma =  conv(Ma, art_kernel);
Ma = Ma(1:Npts);

% label input to the tissue (after the bolus arrival time)
inMa = conv(inMa, art_kernel2);
inMa = inMa(1:Npts);

% Z magnetization in the arteries
Ma =  mtis0*(1 - 2*alpha*Ma);         % stays in the arteries
inMa = mtis0*(1 - 2*alpha*inMa);     % goes into the tissue

aq = 1;

timevec = linspace(0, duration, Npts);

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
    plot(timevec,Ma,'r')
    plot(timevec,inMa,'g')
    hold off
    axis tight
    
    subplot(313)
    plot(obs)
    
end

if isempty(data)
    result = obs;  % Just the signal
else
    data = data(:);
    result = norm(obs(:)- data(:));
    result = result^2;
    %result = obs * conj(data);  % the inverse of the inner product is the cost function
    %result = 1/result;
    LB = [0.0001 ;  0; 0.2; 0.2; 0.35; deg2rad(20)];
    UB = [0.03 ; 0.05; 3; 3; 2; deg2rad(70)];
    UBviolation = find(UB - parms_vector < 0);
    if ~isempty(UBviolation)
        result = result *10;
        %fprintf('\nUB violation');
    end
    LBviolation = find(LB - parms_vector > 0);
    if ~isempty(LBviolation)
        result = result *10;
        %fprintf('\nLB violation');
    end
end
return
