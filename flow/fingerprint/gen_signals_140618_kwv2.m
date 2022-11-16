function obs = gen_signals_140618_kwv2(parms, timing_parms, dofigs,doSub)
% function obs = gen_signals_140328(parms, timing_parms, dofigs,doSub)
% set up parms for pulse sequence and known constants
% ... time units are in seconds

doRK = 1;  % do Runge-kutta method
dt = 0.01;
r1a = 1/1.67;

% get parms from the input structure
mtis0 = parms.mtis0;
f = parms.f;
cbva = parms.cbva;
transit = parms.transit;
kfor = parms.kfor;
r1tis = parms.r1tis;
beta = parms.beta;
Disp = parms.Disp;

if isstruct(timing_parms)
    
    PID =       timing_parms.PID;
    label_dur = timing_parms.label_dur;
    TR =        timing_parms.TR;
    imTR =      timing_parms.imTR;
    TE =        timing_parms.TE;
        
end

N = length(label_dur); % this defines the number of images:
label_start = cumsum(TR);
label_start = [0 label_start(1:end-1)];
aqtimes = cumsum(TR) - (imTR - TE) ;  
duration = sum(TR)+2;

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

%%%%%%%%%%%%%%%%
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
% ttmp(ttmp<0) = 0;
% art_kernel = (Disp^(1+Disp*transit)/gamma(1+Disp*transit)).* ttmp.^(transit*Disp) .* exp(-Disp *ttmp);
% art_kernel = (ttmp).^(transit * Disp) .* exp(-Disp *(ttmp));
%
% art_kernel(ttmp < transit)=0;
% art_kernel = abs(art_kernel);
% art_kernel = art_kernel / ((sum(abs(art_kernel))));  %normalize the mass in the dispersion
% art_kernel = art_kernel .* decay;  


%-------------%

inp = zeros(size(labelfun));
% inp(labelfun==1) = 1;  % in this case we get arterial spin saturation

for n=1:N
    ind1a = floor((label_start(n)+ transit)/dt + 1);
    ind2a = floor((label_start(n) + label_dur(n) + transit)/dt);
    if mod(n,2)==0
        inp(ind1a:ind2a) = 1;
    else
        inp(ind1a:ind2a) = 0;
    end
end
% inp = conv(inp, art_kernel);

% inp is the concentration of label in the artery
% calculate the magnetization in the arterial compartment
% assuming that it is 100% of the voxel and
% the inversion efficiency was 90%
ma = mtis0 *(1 - 2*0.90* inp * exp(-transit*r1a));

% remove padding from both ends
%ma = ma(length(ttmp):end);
ma = ma(1:length(labelfun));
%ma = ma * mtis0/max(ma);


inma = ma;

% now the fun part: calculate the tissue function
% make a magnetization transfer function

mtr = kfor * ~~(labelfun);  % we always get the same RF power and MT whether there is label or not.
m = mtis0 .* ones(size(labelfun));

for n=2:length(t)
    %     if n==3 %Tissue Inversion
    %         m(n-1)=-m(n-1);
    %     end
    
    
    % the modified Bloch equation has
    % t1 decay,  magnetization transfer, inflow, outflow
    
    if ~doRK
        % Euler method:
        %
        dm = (mtis0 - m(n-1))*r1tis - mtr(n-1)*m(n-1) + f*inma(n-1) - f*m(n-1)/0.9;
        m(n) = m(n-1) + dm*dt;
        %
    else
        % Runge-Kutta method:
        %
        tmp = m(n-1);
        k1 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/0.9;
        
        tmp = m(n-1) + k1/2;
        k2 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/0.9;
        
        tmp = m(n-1) + k2/2;
        k3 = (mtis0 - tmp)*r1tis - mtr(n-1)*tmp + f * inma(n-1) - f*tmp/0.9;
        
        tmp = m(n-1) + k3;
        k4 = (mtis0 - tmp)*r1tis - mtr(n)*tmp + f * inma(n) - f*tmp/0.9;
        
        m(n) = m(n-1) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        %
    end
    
    % the Mz gets reset by the sampling RF pulses
    if satfun(n-1)==1
        m(n) =  m(n-1)*cos(beta);
        %ma(n) = 0;
    end
    
    
end

% adjust ma to reflect the blood volume, not just for a 100% voxel
ma = cbva*ma;

obs = m(aqs) + ma(aqs);%%%%%%%%%%%%%%%%%%%%%%%

% Use this if there are DUMMY SCANS
% obs(1) = obs(3);


% modification:  signals are complex.
% Arterial and tissue compartments will acquire different amounts of
% phase during readout


m_phi = 0 ; % 0.1;
ma_phi = 0; % 2;

sigs = ma(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
asub = con-tag;

sigs = m(aqs);
con = sigs(1:2:end);
tag = sigs(2:2:end);
tsub = con-tag;

tot_sub = asub + tsub;

% Let's works with the subtractions
if doSub==1
    obs = tot_sub;
end
obs = obs * sin(beta);

if dofigs
%     figure(10)
%     plot(art_kernel)
    
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
    
    plot(t,ma,'r')
    hold on
    plot(t(aqs),ma(aqs),'g*')
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