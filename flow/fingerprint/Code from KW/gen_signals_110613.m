function obs = gen_signals_110613(parms, getAQfromfile, dofigs)
% set up parms for pulse sequence and known constants
% ... time units are in seconds

if nargin==0
    f =         0.01;
    mtis0 =     1;%0e5 ;
    cbva =      0.05 ;
    transit=    1 ;
    kfor =      1e-4 ;
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

% label_dur=sort(repmat([400:50:1000 1500:500:3000]/1000,[1,16]));
% label_dur=sort(repmat([1:0.1:2.5 3.5],[1,16]));
% label_dur=sort(repmat([0.5:0.1:2.0 3],[1,16]));
% label_dur=sort(repmat([1:0.1:2 3],[1,16]));
% label_dur=sort(repmat([600:100:2100 3300]/1000,[1,16]));
% tr = 1000 + (75*3)* [0:11]
% t_tag = tr - 800
% t_delay = 10

% TR = sort(repmat(((1000 + (100*3)* [0:11])/1000),[1 16])) ; %
% 
% label_dur = TR - 800/1000 ; 
TR = sort(repmat((1200 + 300*[0:11])./1000,[1 16]));
label_dur = TR - 1100/1000;
% t_delay = 10ms 

 


N=length(label_dur);

% TR = label_dur + (600/1000);
duration = sum(TR);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% in these Aq periods, the label will saturate spins, not
% invert.
% mt_aqs = randperm(N);
% mt_aqs = mt_aqs(1:end/4);


aqfun = zeros(N,1);

% labeling function can have 3 values: off(0), label (1), control(-1)
t = 0:dt:duration;


labelfun = zeros(length(t),1);
for n=1:N
    if n==1 || n==2
        if mod(n,2)==0
            label_start(n)=TR(1);
            ind1 = round(label_start(n)/dt);
            ind2 = round((label_start(n) + label_dur(n))/dt);
            labelfun(ind1:ind2) = -1;
        else
            control_start(n)=0;
            ind1 = 1;
            ind2 = round((control_start(n) + label_dur(n))/dt);
            labelfun(ind1:ind2) = 1;
        end
    else
        if mod(n,2)==0
            label_start(n)=control_start(n-1)+TR(n-1);
            ind1 = round(label_start(n)/dt);
            ind2 = round((label_start(n) + label_dur(n))/dt);
            labelfun(ind1:ind2) = -1;
        else
            control_start(n)=label_start(n-1)+TR(n-1);
            ind1 = round(control_start(n)/dt);
            ind2 = round((control_start(n) + label_dur(n))/dt);
            labelfun(ind1:ind2) = 1;
        end
    end
end

aqtime(1)= label_dur(1) + 2*dt + (10/1000);
for n=2:N
    aqtime(n)=aqtime(n-1)+ TR(n);
end
aqfun=zeros(size(t));
aqfun(round(aqtime/dt)) = 1;
aqs = find(aqfun);

% aqfun2=ones(size(t));
% aqfun2(round(aqtime/dt))=0;
% figure;plot(t(1:17741),labelfun(1:17741),'k-',t(aqs),aqfun2(aqs),'r*');ylabel('L(t)');xlabel('Time (s)')
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

for n=2:length(t)-1
    % the modified Blocj equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f*ma(n-1) - f*m(n-1)/0.9;

    m(n) = m(n-1) + dm*dt;

    if aqfun(n-1)==1
        m(n) = cos(deg2rad(90)).*m(n);
        ma(n) = cos(deg2rad(90)).*ma(n);
    end

end


% adjust ma to reflect the blood volume, not just for a 100% voxel

ma = cbva*0.9*ma';

%m = m + ma;
if f>0 && transit>0
    obs = m(aqs) + ma(aqs);%%%%%%%%%%%%%%%%%%%%%%%
else
    obs=m(aqs);
    ma(aqs)=0;
end


% modification:  signals are complex.   Assume flip angle beta is pi/2
% Arterial and tissue compartments will acquire different amounts of
% phase during readout
beta=deg2rad(20);
m_phi = 0.1;
ma_phi = 2;
obs = m(aqs)*sin(beta)*exp(-i*m_phi) + ma(aqs)*sin(beta)*exp(-i*ma_phi);
obs = abs(obs(3:end));
% obs=abs(obs./min(obs(:)));
%obs = ma(aqs);

% figure;plot(t(1:17741),labelfun(1:17741),'k-',t(aqs),aqfun2(aqs),'r*');ylabel('L(t)');xlabel('Time (s)')
% 
% figure;subplot(4,1,1),plot(t(1:4000),labelfun(1:4000));ylabel('L(t)');xlabel('Time (s)')
% subplot(4,1,2),plot(t(1:4000),ma(1:4000),'-r',t(aqs(1:find(aqs>4000,1)-1)),ma(aqs(1:find(aqs>4000,1)-1)),'xr');ylabel('M_a(t)');xlabel('Time (s)')
% subplot(4,1,3),plot(t(1:4000),m(1:4000),'-k',t(aqs(1:find(aqs>4000,1)-1)),m(aqs(1:find(aqs>4000,1)-1)),'xk');ylabel('M_t_i_s(t)');xlabel('Time (s)')
% subplot(4,1,4),plot(t(aqs(1:find(aqs>4000,1)-1)),obs(1:find(aqs>4000,1)-1));ylabel('M_o_b_s(t)');xlabel('Time (s)')

if dofigs
    figure
    subplot(411)
    area(t,labelfun)
    hold on
    plot(t,tmp_inp)
%     plot(t,bsfun,'r')
%     plot(t(aqs),bsfun(aqs),'g*')
    title('label/control/BS pulses')
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


    figure
    subplot(211)
    plot(abs(obs)); title('signal magnitude')
    subplot(212)
    plot(angle(obs)); title('signal phase')

    figure
    plot(abs(obs))
   
    drawnow
end

% sigs = ma(aqs);
% con = sigs(1:2:end);
% tag = sigs(2:2:end);
% asub = con-tag;
% 
% sigs = m(aqs);
% con = sigs(1:2:end);
% tag = sigs(2:2:end);
% tsub = con-tag;
% 
% figure
% plot( asub); hold on; plot(tsub,'r')

end