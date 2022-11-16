function obs = gen_signals_140127(parms, timing_parms, dofigs)
% set up parms for pulse sequence and known constants
% ... time units are in seconds
% this versioni is intended for testing at T1 fit to phantom data

if nargin==0

    
    f =         0.015;
    mtis0 =     10e5 ;
    cbva =      0.02 ;
    transit=    1.2 ;
    kfor =      0.2; % 1e-2 ;
    r1tis =     1/1.4  ;
    beta =      15*pi/180 ; % flip angle in radians
    L = 1;
    %
    
%     f = 0.01 
%     cbva =0.0323 
%     transit = 1.2
%     kfor = 0.2632 
%     r1tis = 1.3971
%     mtis0 =     10e5 ;
%     beta =      10*pi/180 ; % flip angle in radians
%     L = 1;
%     %

    
    parms = struct( ...
        'mtis0', mtis0,...
        'f', f , ...  % perfusion in ml/s/g
        'cbva' ,  cbva , ...
        'transit', transit,...
        'kfor', kfor, ...
        'r1tis', r1tis, ...
        'beta', beta);
    
    getAQfromfile=0;
    dofigs = 1;
end

dt = 0.01;
Disp = 0.1/dt;
Disp = 2;  % more realistic
r1a = 1/1.6;

% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;
beta = parms.beta;

mt_aqs = [];

PID = 0.01;
TR = sort(repmat( [1:0.15:3.85], [1,4]));
label_dur = sort(repmat([0.1 : 0.15 : 2.95 ], [1,4])); %  ]TR - 0.600 - 0.300;

%% reverse the order of the TRs, etc 
%{
oldorder = 1:length(TR);
neworder = oldorder(end:-1:1);

label_dur = label_dur(neworder);
TR = TR(neworder);
%}
%%

TR =        sort(repmat( linspace(2.2, 4, 40) , [1,4]));
label_dur = sort(repmat( linspace(0.2, 1.5, 40) , [1,4])); %  ]TR - 0.600 - 0.300;
PID =       sort(repmat( linspace(0.2, 1.8, 40) , [1,4]));

PID(3:8:end) = 1.6;
PID(4:8:end) = 1.6;

label_dur(1:8:end) = 0;
label_dur(2:8:end) = 0;

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

aqtimes = cumsum(TR) - 0.6 ;  % acquisition happens 10 ms after labeling
duration = sum(TR)+2;


% this defines the number of images:
N = length(label_dur);

label_start = cumsum(TR) - label_dur - PID - 0.6;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% labeling function can have 3 values: off(0), label (1), control(-1)
t = 0:dt:duration;


labelfun = zeros(length(t),1);
aqfun = zeros(size(t));
satfun = zeros(size(t));


for n=1:N
    % begin label
    ind1 = floor(label_start(n)/dt + 1);
    %end label
    ind2 = floor((label_start(n) + label_dur(n))/dt);
    
    if mod(n,2)==0
        labelfun(ind1:ind2) = -1;
    else
        labelfun(ind1:ind2) = 1;
    end
    
    % acquire image 100 ms after end of label
    ind3 = floor(aqtimes(n)/dt);
    aqfun(ind3) = 1;

    % Magnetization is tipped during whole acquisition period
    % saturation function to describe when the magnetization is tipped for imaging 
    ind4 = ind3 + [0:0.03:0.6]/dt;
    satfun(floor(ind4)) = 1;

    
end
aqfun = aqfun(1:length(labelfun));

aqs = find(aqfun);


%%% remove background suppression  ****

%%%%%%%%%%%%%%%%
% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,5, 5/dt)';
art_kernel = (ttmp-transit).* exp(-Disp *(ttmp-transit));
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion

decay = ones(size(ttmp));
decay = exp(-ttmp*r1a);
art_kernel = art_kernel .* decay;



inp = zeros(size(labelfun));
inp(labelfun==1) = 1;  % in this case we get arterial spin saturation
%inp(labelfun==0.5) = 0;  % in this case we get arterial spin saturation

% convolve with dipersion kernel
% in order to obtain the arterial signal
% assume that the arterial volume is 100% for input function
ma =  mtis0 *(1 - 2*0.8* conv(inp, art_kernel));

% remove padding from both ends
%ma = ma(length(ttmp):end);
ma = ma(1:length(labelfun));
%ma = ma * mtis0/max(ma);

% there is a short lag between arriving at the voxel and getting into the
% tissue.  It's from 200 to 500 ms.
trans2 = 0.200 ;
inma = zeros(size(ma));
inma(trans2/dt +1: end) = ma(1: end - trans2/dt) ;

% now the fun part: calculate the tissue function
% make a matgnetization transfer function

mtr = kfor * ~~(labelfun);  % we always get the same RF power and MT whether there is label or not.
m = mtis0 * ones(size(labelfun));

for n=3:length(t)
    % the modified Blocj equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f * inma(n-1) - f*m(n-1)/0.9;

    m(n) = m(n-1) + dm*dt;
    
    % the Mz gets reset by the sampling RF pulses
    if satfun(n-2)==1
        m(n) =  0 ; % m(n-1)*cos(beta);
        %ma(n) = 0;
    end
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
ma = cbva*0.9*ma;

obs = m(aqs) + ma(aqs);%%%%%%%%%%%%%%%%%%%%%%%

% Use this if there are DUMMY SCANS
obs(1) = obs(3);


% modification:  signals are complex.
% Arterial and tissue compartments will acquire different amounts of
% phase during readout


m_phi = 0 ; % 0.1;
ma_phi = 0; % 2;

% obs = m(aqs)*sin(beta)*exp(-i*m_phi) + ma(aqs)*sin(beta)*exp(-i*ma_phi);
%obs = obs(3:end);
%obs = ma(aqs);

obs = obs * sin(beta);

% Let's works with the subtractions
obs = obs(2:2:end) - obs(1:2:end);


sigs = ma(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
asub = con-tag;

sigs = m(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
tsub = con-tag;

tot_sub = asub + tsub;

if dofigs
    figure(10)
    plot(art_kernel)

    figure(1)
    subplot(411)
    area(t,labelfun)
    hold on
    %    plot(t,tmp_inp)
    %     plot(t,bsfun,'r')
    %     plot(t(aqs),bsfun(aqs),'g*')
    title('label/control/BS pulses')
    hold on
    plot(t(aqs),labelfun(aqs),'g*')
    hold off
    
    subplot(412)
    
    plot(t,ma,'r')
    hold on
    plot(t(aqs),ma(aqs),'g*')
    hold off
    title('arterial magnetization')
    
    subplot(413)
    plot(t,m)
    hold on
    plot(t(aqs),m(aqs),'g*')
    hold off
    title('tissue magnetization')
    
    subplot(414)
    plot(angle(obs))
    
    title('phase of aquired signal')
    
    
    figure(5)
    subplot(211)
    plot(abs(obs)); title('signal magnitude')
    subplot(212)
    plot(angle(obs)); title('signal phase')
    
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
