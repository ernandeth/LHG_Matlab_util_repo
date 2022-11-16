function obs = gen_signals(parms, getAQfromfile, dofigs)
% set up parms for pulse sequence and known constants
% ... time units are in seconds

if nargin==0
    f =         0.01;
    mtis0 =     10e2 ;
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
    close all
end

duration = 90;
dt = 0.01;
Disp = 0.1/dt;
r1a = 1/1.6;
TR = 4;


% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;

doBS = 1;

if ~getAQfromfile
%     fprintf('\nGenerating a new random acquisition sequence ...');
    option=5;
    
    mt_aqs = [];
    
    switch(option)
            case 5
            % a sequential increase in labeling times with gap (PLD) in between
            % followed by the same thing without PLD

            label_dur = 0:0.3:2.5;

            % control and tag have the same duration and are acquired
            % alteranting
            % starttime is the start of the sequence
            starttime = zeros(1, length(label_dur)*2);
            starttime(3:2:end) = label_dur(2:end) + 0.2;
            starttime(2:2:end) = label_dur + 0.2;
            
            % with PLD
            aqtime = cumsum(starttime + 1.5);
            % without the PLS
            aqtime2 = aqtime(end) + label_dur(end) + 1.5 + cumsum(starttime + 0.2);
            
            % merge the two sequences into one longer one.
            aqtime = [aqtime aqtime2];
            label_dur = [label_dur label_dur];

            
            bstime0 = zeros(size(aqtime));
            bstime1 = 2.1 * ones(size(aqtime));
            bstime2 = 3.0 * ones(size(aqtime));
            
        
            
    end
      
      
    save AQparms label_dur  bstime1 bstime2 bstime0 aqtime mt_aqs TR
else
    fprintf('\rReading acquisition sequence from file...');
    load AQparms
end

labelfun = zeros(duration/dt,1);
bsfun = zeros(duration/dt,1);
aqfun = zeros(duration/dt,1);
t = linspace(0,duration,length(labelfun));

% labeling function can have 3 values: off(0), label (1), control(-1)
%

for n=1:length(label_dur)-1    
    
    ind1 = aqtime(2*n-1) ;
    ind2 = aqtime(2*n-1) + label_dur(n);

    ind3 = aqtime(2*n-1) + bstime0(n);
    ind4 = aqtime(2*n-1) + bstime1(n);
    ind5 = aqtime(2*n-1) + bstime2(n);
    
   
    labelfun( round(ind1/dt): round(ind2/dt) ) = 1;
    bsfun(round([ind3 ind4 ind5]/dt)) = 1;
    
    
    ind1 = aqtime(2*n) ;
    ind2 = aqtime(2*n)  + label_dur(n);
    
    ind3 = aqtime(2*n) + bstime0(n);
    ind4 = aqtime(2*n) + bstime1(n);
    ind5 = aqtime(2*n) + bstime2(n);
    
    labelfun( round(ind1/dt): round(ind2/dt) ) = -1;
    bsfun(round([ind3 ind4 ind5]/dt)) = 1;

    
end


bsfun(2) = 1;
aqfun(round(aqtime/dt)) = 1;
aqs = find(aqfun);


% interesting Background suppression scheme:
bsfun(:) = 0;
inds = round([0.1:0.7:duration]/dt);
bsfun(inds) = 1;


%%% remove background suppression  ****
if ~doBS
    bsfun(:) = 0;
end
%%%%%%%%%%%%%%%%

% form an arterial dispersion kernel (only need to do this once)
ttmp = linspace(0,5, 5/dt)';
art_kernel = (ttmp-transit).* exp(-Disp *(ttmp-transit));
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);


inp = ones(size(labelfun));
inp(labelfun==1) = 1 - 2*0.8* exp(-transit*r1a);
%inp(labelfun==0.5) = 1 - exp(-transit*r1a);  % in this case we get arterial spin saturation
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

%kfor =0;

mtr = kfor * ~~(labelfun); 
m = mtis0 * ones(size(t));

for n=2:length(t)
    
    dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f*ma(n-1) - f*m(n-1);
    
    m(n) = m(n-1) + dm*dt;
    
    if bsfun(n)==1
        m(n) = -m(n);
        ma(n) = -ma(n);
    end
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
ma = cbva*0.9*ma'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coefficient place?

obs = m(aqs) + ma(aqs);

beta = pi/2;
m_phi = 0.1;
ma_phi = 2;
% obs = m(aqs)*sin(beta)*exp(-i*m_phi) + ma(aqs)*sin(beta)*exp(-i*ma_phi);
obs = m(aqs) + ma(aqs);

obs = obs(3:end);


if dofigs
    figure(1)
    subplot(311)
    area(t,labelfun)
    hold on
    plot(t,tmp_inp)
    plot(t,bsfun,'r')
    plot(t(aqs),bsfun(aqs),'g*')
    title('label/control/BS pulses')
    hold off
    
    subplot(312)

    hold on
    plot(t,ma,'r')
    plot(t(aqs),ma(aqs),'g*')
    hold off
    title('arterial magnetization')
    
    subplot(313)
    plot(t,m,'b')
    hold on
    plot(t(aqs),m(aqs),'g*')
    hold off
    title('tissue magnetization')
    
%    
%     figure(5)
%     subplot(211)
%     plot(abs(obs)); title('signal magnitude')
%     subplot(212)
%     plot(angle(obs)); title('signal phase')
%     
    
    drawnow
end



end 
