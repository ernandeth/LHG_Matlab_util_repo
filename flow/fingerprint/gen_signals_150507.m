function obs = gen_signals_150507(parms, timing_parms, dofigs,doSub)
% function obs = gen_signals_150507(parms, timing_parms, dofigs,doSub)
%
% set up parms for pulse sequence and known constants
% ... time units are in seconds
% assumes that any timing errors in the scannes are already accounted for
% in teh calculation of the timnig_parms
% 
% this version groups the images in 3 tags in a row followed by 3 controls
% in a row ... etc
% 
%
%%
if nargin==0
    % for testing purposes, here are some default parameters:
    f =         0.01;
    f =         60 /6000;
    mtis0 =     1 ;
    cbva =      0.02 ;
    transit =   1.2 ;
    kfor =      1e-2 ;
    r1tis =     1/1.4  ;
    beta =      90*pi/180 ; % flip angle in radians
    Disp =      30;
    L = 1;
    %
    
    doSub = 0;
    
    parms = struct( ...
        'mtis0', mtis0,...
        'f', f , ...  % perfusion in ml/s/g
        'cbva' ,  cbva , ...
        'transit', transit,...
        'kfor', kfor, ...
        'r1tis', r1tis, ...
        'beta', beta, ...
        'Disp', Disp);
    
    getAQfromfile=0;
    dofigs = 1;
    
    
%         TR =        sort(repmat( linspace(0.2, 5, 40) , [1,2]));
%         label_dur = sort(repmat( linspace(0.2, 4, 40) , [1,2])); %  ]TR - 0.600 - 0.300;
%         PID =       sort(repmat( linspace(0.2, 0.8, 40) , [1,2]));
%         t_adjust = 0.1*ones(size( PID )); 
%         
    
    
%     n=1;
%     for myta = [0.5 1]
%         for mypid = [1.6 0.2 1.6]
%             for myld = sort(repmat( linspace(0.2, 2, 10) , [1,2]));
%                 t_adjust(n) = myta;
%                 
%                 t_adjust(n) = 4 - myld - mypid - 0.03;
%                 
%                 label_dur(n) = myld;
%                 PID(n) = mypid;
%                 n = n+1;
%             end
%         end
%     end
    
%    TR = label_dur + PID + 0.03 + t_adjust;
    
    % timing for PCASL :
    PID  = 1.6 * ones(10,1);
    label_dur = 2*ones(10,1);
    t_adjust = 0.4*ones(10,1);
%     
    timing_parms.PID = PID';
    timing_parms.label_dur = label_dur';
    timing_parms.t_adjust = t_adjust';
    timing_parms.order = 1;  % control-tag
    
    timing_parms.Nlabel_group = 1;  % Number of tags in a row  (added this on 5/7/15)
    tmp = ones(size(t_adjust));
    tmp(2:2:end) = -1;
    timing_parms.isLabel = tmp';
    
end
%%

doRK = 0;  % do Runge-kutta method
dt = 0.01; % seconds

% Disp = 0.1/dt;
% Disp = 10;  % more realistic
r1a = 1/1.67;
lambda = 0.9;

% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;
beta = parms.beta;
Disp = parms.Disp;

mt_aqs = [];



if isstruct(timing_parms)
    
    PID =       timing_parms.PID;
    label_dur = timing_parms.label_dur;
    t_adjust =        timing_parms.t_adjust;
    Nlabel_group =    timing_parms.Nlabel_group;
    order =     timing_parms.order;
    isLabel = timing_parms.isLabel ;
end

aqwindow = 0.032 ;

% acquisition takes 30 ms. from RF to end of crusher
TR = t_adjust + PID  + label_dur + aqwindow;

aqtimes = cumsum(TR) - aqwindow  ;
duration = sum(TR)+2;

% this defines the number of images:
N = length(label_dur);

label_start = aqtimes - label_dur - PID ;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% labeling function can have 3 values: off(0), label (1), control(-1)
t = 0:dt:duration;


labelfun = zeros(length(t),1);
aqfun = zeros(size(t));
satfun = zeros(size(t));

ctr = 1;

for n=1:N
    
   
    % begin label
    ind1 = round((label_start(n))/dt ) +1;
    %end label
    ind2 = round((label_start(n) + label_dur(n))/dt)+1;
    
    % The default order is tag-control
    % if mod(n,2)==0          % even->control
    %    labelfun(ind1:ind2) = -1;
    % else                            % odd -> tag
    %    labelfun(ind1:ind2) = 1;
    % end
             
    
     labelfun(ind1:ind2) = isLabel(n); % tag
     
     % if ctr <= Nlabel_group ,         labelfun(ind1:ind2) = -1;     end % control
     % if ctr == 2*Nlabel_group ,     ctr=0;  end
         
     ctr = ctr + 1;
     
    % acquire image
    ind3 = round(aqtimes(n)/dt);
    aqfun(ind3) = 1;
    
    % Magnetization is tipped during whole acquisition period
    % saturation function to describe when the magnetization is tipped for imaging
    %  this is for multi-slice 3D imaging acquisition
    ind4 = ind3 ; %   + [0:0.001:0.02]/dt;
    satfun((ind4)) = 1;
    
    
end
aqfun = aqfun(1:length(labelfun));

aqs = find(aqfun);

if timing_parms.order==2  
    labelfun = -labelfun;
end


%%% remove background suppression  ****

%%%%%%%%%%%%%%%%
% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,10, 5/dt)';
decay = ones(size(ttmp));
decay = exp(-ttmp*r1a);

%{
art_kernel = (ttmp-transit).* exp(-Disp *(ttmp-transit));
%art_kernel = (ttmp).^(transit * Disp) .* exp(-Disp *(ttmp));
%
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
art_kernel = art_kernel .* decay;

%}
%--------- trying out new kernel from Chappell's paper


%ttmp = ttmp - transit;
% ttmp(ttmp<0) = 0;
art_kernel = ttmp.^(transit*Disp) .* exp(-Disp *ttmp);

art_kernel(ttmp <= 0)=0;
art_kernel = abs(art_kernel);
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
art_kernel = art_kernel .* decay;

%-------------%

inp = zeros(size(labelfun));
inp(labelfun==1) = 1;  % in this case we get arterial spin saturation
%inp(labelfun==0.5) = 0;  % in this case we get arterial spin saturation

%
% ignore dispersion - consider only T1 decay during transit:
ma_tmp = mtis0 *(1 - 2*0.90* inp * exp(-transit * r1a) );
ma_tmp =  circshift(ma_tmp, round(transit/dt));
ma_tmp = ma_tmp(1:length(labelfun));

%ma = ma_tmp;

% convolve with dipersion kernel
% in order to obtain the arterial concentration of the label
% assume that the arterial volume is 100% for input function
inp = conv(inp, art_kernel);


% inp is the concentration of label in the artery
% calculate the magnetization in the arterial compartment
% assuming that it is 100% of the voxel and
% the inversion efficiency was 90%

ma =  mtis0 *(1 - 2*0.90* inp);

% .... OR .....
% Neglect dispersion effects:  just use the shifted/decayed input function
%ma = ma_tmp;
%

% remove padding from both ends
% ma = ma(length(ttmp):end);
ma = ma(1:length(labelfun));


% there is a short lag between arriving at the voxel and getting into the
% tissue.  It's from 200 to 500 ms.
% trans2 = 0.00 ;
% inma = zeros(size(ma));
% inma(trans2/dt +1: end) = ma(1: end - trans2/dt) ;
% ... never mind ....
inma = ma;

% now the fun part: calculate the tissue function
% make a magnetization transfer function

mtr = kfor * ~~(labelfun);  % we always get the same RF power and MT whether there is label or not.
m = mtis0 * ones(size(labelfun));

%mtr(:) = 0;

% Using Bloch equations to calculate magnetization in tissue, takes into
% account MR sampling, T1 , magnetization transfer ...etc.
for n=2:length(t)
    
    
    % the modified Bloch equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    
    
    if doRK
        % Runge-Kutta method:
        %
        tmp = m(n-1);
        k1 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/lambda;
        
        tmp = m(n-1) + k1/2;
        k2 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/lambda;
        
        tmp = m(n-1) + k2/2;
        k3 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/lambda;
        
        tmp = m(n-1) + k3;
        k4 = (mtis0 - tmp)*r1tis - mtr(n)*tmp + f * inma(n) - f*tmp/lambda;
        
        m(n) = m(n-1) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        %
        
    else
        % Euler method:
        %
        dm = (mtis0 - m(n-1))*r1tis  - mtr(n-1)*m(n-1) + f * inma(n-1) - f*m(n-1)/lambda;
        
        %dm = (mtis0 - m(n-1))*r1tis; % Just T1 !!!
        
        m(n) = m(n-1) + dm*dt;
        
    end
    
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if aqfun(n-1)==1
        m(n) =  m(n-1)  * cos(beta);
        %ma(n) = ma(n-1)*cos(beta);
    end


    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
m = (1-cbva)*m;
ma = cbva*ma;
obs = m(aqs) + ma(aqs);
obs = obs * sin(beta);

obs = m(aqs) * sin(beta);


% Use this if there are DUMMY SCANS
% obs(1) = obs(3);


sigs = ma(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
asub = con-tag;

sigs = m(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
tsub = con-tag;



if doSub==1
    tot_sub = obs(1:2:end) -  obs(2:2:end);
    
    obs = tot_sub;
end

% Sanity check:  Using convolution to calculate concentration of label in tissue:
%{
ttmp = linspace(0,10, 10/dt)';
ret = exp( -(f/lambda + r1tis) * ttmp);
mt =  2*0.9*mtis0 * f * conv(inp, ret ) * dt;

mt = mt(1:length(t));
inp = inp(1:length(t));
figure(13)
plot(t,mt); hold on
plot(t,inp,'r');
plot(t(aqs),mt(aqs),'g*')
hold off
sub = diff(mt(aqs));
signal = mean(abs(sub))

%}

if dofigs
    figure(10)
    plot(ttmp, art_kernel)
    
    figure(1)
    subplot(311)
    area(t,labelfun)
    hold on
    
    %    plot(t,tmp_inp)
    %     plot(t,bsfun,'r')
    %     plot(t(aqs),bsfun(aqs),'g*')
    
    title('label/control/BS pulses')
    hold on
    plot(t(aqs),labelfun(aqs),'g*')
    axis tight
    hold off
    
    subplot(312)
    
    plot(t,ma/cbva,'r')
    hold on
    plot(t,ma_tmp','b')
    plot(t(aqs),ma(aqs)/cbva,'g*')
    axis tight
    hold off
    title('arterial magnetization')
    
    subplot(313)
    plot(t,m)
    hold on
    plot(t(aqs),m(aqs),'g*')
    hold off
    title('tissue magnetization')
    axis tight
    
    %     subplot(414)
    %     plot(angle(obs))
    %
    %     title('phase of aquired signal')
    
    
    %     figure(5)
    %     subplot(211)
    %     plot(abs(obs)); title('signal magnitude')
    %     subplot(212)
    %     plot(angle(obs)); title('signal phase')
    %
    
    figure (6)
    subplot(211)
    plot(ma(aqs),'r');
    subplot(212)
    plot(ma(aqs) + m(aqs),'k')
    hold on
    plot(m(aqs),'b')
    legend('Both', 'Tissue Only')
    grid on
    hold off
    
    drawnow
end

end
