function obs = gen_signals_140707(parms, timing_parms, dofigs,doSub)
% function obs = gen_signals_140627(parms, timing_parms, dofigs,doSub)
% set up parms for pulse sequence and known constants
% ... time units are in seconds
% this versioni is intended for testing at T1 fit to phantom data

if nargin==0
    
    f =         0.0133;
    mtis0 =     1 ;
    cbva =      0.02 ;
    transit =   1.2 ;
    kfor =      0.2; % 1e-2 ;
    r1tis =     1/1.4  ;
    beta =      90*pi/180 ; % flip angle in radians
    Disp =      30;
    L = 1;
    %
    
    doSub = 1;
    
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
    psdseqtime = 0.0346;

    
    TR =        sort(repmat( linspace(0.2, 5, 40) , [1,2]));
    label_dur = sort(repmat( linspace(0.2, 4, 40) , [1,2])); %  ]TR - 0.600 - 0.300;
    PID =       sort(repmat( linspace(0.2, 0.8, 40) , [1,2]));
%     label_dur = load('t_tags.txt');
%     PID = load('t_delays.txt');
%     TR = label_dur + PID + psdseqtime + 0.035;
    
    %PID(3:8:end) = 1.6;
    %PID(4:8:end) = 1.6;
    
    %label_dur(1:8:end) = 0;
    %label_dur(2:8:end) = 0;
    
    % timing for a turbo curve:
    % TR = label_dur + 0.2;
    % PID (:) = 0.01;
    
    % timing for PCASL :
    TR(:) = 5.;
    PID (:) = 1.6;
    label_dur(:) = 3;
    
    timing_parms.PID = PID;
    timing_parms.label_dur = label_dur;
    timing_parms.TR = TR;
end

doRK = 0;  % do Runge-kutta method
dt = 0.01;
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

% PID = 0.01;
% TR = sort(repmat( [1:0.15:3.85], [1,4]));
% label_dur = sort(repmat([0.1 : 0.15 : 2.95 ], [1,4])); %  ]TR - 0.600 - 0.300;

%% reverse the order of the TRs, etc
%{
oldorder = 1:length(TR);
neworder = oldorder(end:-1:1);

label_dur = label_dur(neworder);
TR = TR(neworder);
%}
%%


%% shuffle the order

% inds2 = randperm(length(TR)); save myorder.dat inds2 -ascii
% inds2 = load('myorder.dat');
%
% PID = PID(inds2);
% label_dur = label_dur(inds2);
% TR = TR(inds2);
%
%%

if isstruct(timing_parms)
    
    PID =       timing_parms.PID;
    label_dur = timing_parms.label_dur;
    TR =        timing_parms.TR;
    
end
aqtimes = cumsum(TR) - 0.07 ;  % acquisition happens 10 ms after labeling
duration = sum(TR)+2;


% this defines the number of images:
N = length(label_dur);

label_start = cumsum(TR)' ; % - label_dur - PID - 0.02;
label_start= [0 label_start(1:end-1)];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% labeling function can have 3 values: off(0), label (1), control(-1)
t = 0:dt:duration;


labelfun = zeros(length(t),1);
aqfun = zeros(size(t));
satfun = zeros(size(t));


for n=1:N
    % begin label
    if n==1
        ind1 = floor(label_start(n)/dt + 1);
    else
        % ind1 = floor((label_start(n) + 0.077)/dt + 1);
        % this gap is already accounted for in the TR
        ind1 = floor((label_start(n))/dt + 1);
    end
    %end label
    ind2 = floor((label_start(n) + label_dur(n))/dt);
    
    if mod(n,2)==0
        labelfun(ind1:ind2) = 1;
    else
        labelfun(ind1:ind2) = -1;
    end
    
    % acquire image 100 ms after end of label
    ind3 = floor(aqtimes(n)/dt);
    aqfun(ind3) = 1;
    
    % Magnetization is tipped during whole acquisition period
    % saturation function to describe when the magnetization is tipped for imaging
    ind4 = ind3 ; %   + [0:0.001:0.02]/dt;
    satfun(floor(ind4)) = 1;
    
    
end
aqfun = aqfun(1:length(labelfun));

aqs = find(aqfun);


%%% remove background suppression  ****

% %%%%%%%%%%%%%%%%
% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,5, 5/dt)';
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


% ttmp = ttmp - transit;
ttmp(ttmp<0) = 0;
art_kernel = ttmp.^(transit*Disp) .* exp(-Disp *ttmp);

art_kernel(ttmp <= 0)=0;
art_kernel = abs(art_kernel);
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion
art_kernel = art_kernel .* decay;

%-------------%

inp = zeros(size(labelfun));
inp(labelfun==1) = 1;  % in this case we get arterial spin saturation
%inp(labelfun==0.5) = 0;  % in this case we get arterial spin saturation
% inp =  circshift(inp, round(transit/dt));
initial_inp=inp;

% convolve with dipersion kernel
% in order to obtain the arterial concentration of the label
% assume that the arterial volume is 100% for input function

inp = conv(inp, art_kernel);
ma =  mtis0 *(1 - 2*0.90* inp);


% inp is the concentration of label in the artery
% calculate the magnetization in the arterial compartment
% assuming that it is 100% of the voxel and
% the inversion efficiency was 90%

%ma =  mtis0 *(1 - 2*0.90* inp*exp(-transit*r1a));

% 
% ignore dispersion - consider only T1 decay during transit:
% ma = mtis0 *(1 - 2*0.90* inp * exp(-transit * r1a) );
% ma =  circshift(ma, round(transit/dt));

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

% mtr(:) = 0;

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
        dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f * inma(n-1) - f*m(n-1)/lambda;       
        
        m(n) = m(n-1) + dm*dt;
    
    end
    
    %  Mz gets tipped toward the xy plane by the sampling RF pulses
    if satfun(n-1)==1
        m(n) =  m(n-1)*cos(beta);
        % ma(n) = 0;
    end
    
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
% ma = cbva*ma;

% obs = m(aqs) + ma(aqs);%%%%%%%%%%%%%%%%%%%%%%%

obs= (1-cbva)*m(aqs) + cbva*ma(aqs);

% Use this if there are DUMMY SCANS
% obs(1) = obs(3);


% modification:  signals are complex.
% Arterial and tissue compartments will acquire different amounts of
% phase during readout

m_phi = 0 ; % 0.1;
ma_phi = 0; % 2;

sigs = cbva*ma(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
asub = con-tag;

sigs =(1-cbva)* m(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
tsub = con-tag;

tot_sub = asub+tsub;

% Let's works with the subtractions
if doSub==1
    obs = tot_sub;
end
obs = obs * sin(beta);

% % Sanity check:  Using convolution to calculate concentration of label in tissue:
% ttmp = linspace(0,10, 10/dt)';
% ret = exp( -(f/lambda + r1tis) * ttmp);
% mt =  2*0.9*mtis0 * f * conv(inp, ret ) * dt;
% 
% mt = mt(1:length(t));
% inp = inp(1:length(t));
% % figure(13)
% % plot(t,mt); hold on
% % plot(t,inp,'r');
% % plot(t(aqs),mt(aqs),'g*')
% hold off
% sub = diff(mt(aqs));
% signal = mean(abs(sub));

%

if dofigs
    figure(10)
    plot(art_kernel)
    
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
    
    plot(t,ma,'b--')
    hold on
    plot(t(aqs),ma(aqs),'g*')
    hold on
    plot(t,initial_inp,'r');
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
    plot(asub,'r'); hold on;
    plot(tsub,'b')
    plot(tot_sub,'k')
    title('Subtractions')
    legend('arterial only','tissue only', 'Both')
    grid on
    hold off
    
    drawnow
end

end
