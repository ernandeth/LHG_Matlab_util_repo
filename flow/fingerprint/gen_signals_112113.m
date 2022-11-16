function obs = gen_signals(parms, getAQfromfile, dofigs)
% set up parms for pulse sequence and known constants
% ... time units are in seconds

if nargin==0
    f =         0.01;
    mtis0 =     10e5 ;
    cbva =      0.05 ;
    transit=    1.2 ;
    kfor =      1e-2 ;
    r1tis =     1/1.4  ;
    L = 1;
    parms = struct( ...
        'mtis0', mtis0,...
        'f', f , ...  % perfusion in ml/s/g
        'cbva' ,  cbva , ...
        'transit', transit,...
        'kfor', kfor, ...
        'r1tis', r1tis);
    getAQfromfile=0;
    dofigs = 1;
end

dt = 0.01;
Disp = 0.1/dt;
r1a = 1/1.6;

% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;

mt_aqs = [];

% select the tagging durations :
label_dur=sort(repmat([1:0.1:2.5 3.5],[1,16]));
TR = label_dur + 0.600;
aqtimes = label_dur + 0.01 ;  % acquisition happens 10 ms after labeling
duration = sum(TR);

% this defines the number of images:
N = length(label_dur);

label_start = cumsum(TR) - TR(1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% labeling function can have 3 values: off(0), label (1), control(-1)
t = 0:dt:duration;

labelfun = zeros(length(t),1);
aqfun = zeros(size(t));

for n=1:N-2
    % begin label
    ind1 = round(label_start(n)/dt+1)
    %end label
    ind2 = round((label_start(n) + label_dur(n))/dt)
    if mod(n,2)==0
        labelfun(ind1:ind2) = -1;
    else
        labelfun(ind1:ind2) = 1;        
    end
    % acquire image
    ind3 = ind2 + 0.01/dt;
    aqfun(ind3) = 1;
    
end

aqs = find(aqfun);


%%% remove background suppression  ****

%%%%%%%%%%%%%%%%
% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,5, 5/dt)';
art_kernel = (ttmp-transit).* exp(-Disp *(ttmp-transit));
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);


inp = ones(size(labelfun));
inp(labelfun==1) = 1 - 2*0.8* exp(-transit*r1a);
inp(labelfun==0.5) = 1 - exp(-transit*r1a);  % in this case we get arterial spin saturation
tmp_inp = inp;

% pad the input function with ones
inp = [ones(size(ttmp)); inp];

% convolve with dipersion kernel
% in order to obtain the arterial signal
% assume that the arterial volume is 100% for input function
ma = mtis0 * conv(inp, art_kernel);

% remove padding from both ends
ma = ma(length(ttmp):end);
ma = ma(1:length(labelfun));
ma = ma * mtis0/max(ma);



% now the fun part: calculate the tissue function
% make a matgnetization transfer function

mtr = kfor * ~~(labelfun);  % we always get the same RF power and MT whether there is label or not.
m = mtis0 * ones(size(t));

for n=7:length(t)
    % the modified Blocj equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f*ma(n-1-5) - f*m(n-1)/0.9;
    m(n) = m(n-1) + dm*dt;
  
    if aqfun(n-1)==1
        m(n) = 0;
        ma(n) = 0;
    end
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
ma = cbva*0.9*ma';
%m = m + ma;

obs = m(aqs) + ma(aqs);%%%%%%%%%%%%%%%%%%%%%%%


% modification:  signals are complex.   Assume flip angle beta is pi/2
% Arterial and tissue compartments will acquire different amounts of
% phase during readout
beta = pi/2;
m_phi = 0.1;
ma_phi = 2;
obs = m(aqs)*sin(beta)*exp(-i*m_phi) + ma(aqs)*sin(beta)*exp(-i*ma_phi);
obs = obs(3:end);
%obs = ma(aqs);

if dofigs
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
    
    hold on
    plot(t,ma,'r')
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
    
    
    drawnow
end

sigs = ma(aqs);     
con = sigs(1:2:end);
tag = sigs(2:2:end);
asub = con-tag;

sigs = m(aqs);     
con = sigs(1:2:end);
tag = sigs(2:2:end);
tsub = con-tag;

figure
plot( asub); hold on; 
plot(tsub,'r')

end
