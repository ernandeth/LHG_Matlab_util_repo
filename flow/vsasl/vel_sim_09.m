%close all
% figure(1)

% simulation step size in milliseconds.
% 1 usec per step. (or 0.001 msec)
dt = 1e-3;
showEvolution = 0;
isControl = 0;
    
T1 = 1700;  %ms
T2 = 165;   %ms  from Qin paper
%        
% T1 = 1400;  %ms
% T2 = 65;   %ms  from Qin paper
  
homogeneity = 'perfect'; % 'perfect'  % 'bad_B0_B1'
refocus_scheme = 'MLEV'; % 'MLEV' % 'DR180'
base_pulse =  'sinc-mod'; %'hard'   % 'sinc'  % 'hard'; % sech, hard, BIR

batch_mode = 0;
exportPulses = 0;

if batch_mode
    homogeneity ='';
end


 
switch homogeneity
    case 'bad_B0'
        off_resonance =    100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1;
        
    case 'bad_B0_B1'
        off_resonance =   100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1.3;
        
    case 'perfect'
        off_resonance = 0%  100e-3 / 42576 ; % (KHz --> Tesla) ... gammabar is in kHz/T
        B1fudge = 1 % 1.1;
end




% default movement parameters
accel = 0 ; % in cm/ms^2
vel = 0  ; % in cm/ms
pos0 = 0 ;  % cm

Mzfinal = [];
Mzfinal_ns = [];


Nsegs = 9;

Gvs =   1.50*1e-4;  % T/cm
%Gvs =   1.00*1e-4;  % T/cm
Gvs =   3*1e-4;  % T/cm

gap = round(0.4/dt);
gap = round(0.2/dt); % in the paper
ramp_len = 0.3/dt ;
%ramp_len = 0.28/dt;
trap = [
    zeros(gap,1);
    [0:ramp_len-1]'/(ramp_len-1);
    [ramp_len-1:-1:0]'/(ramp_len-1);
    zeros(gap,1);
    ]*Gvs;
trap_len = length(trap);

B1_area180 = 0.117e-4; % Tesla*ms


switch base_pulse
    
    case 'sinc'
        
        mySinc_dur = 2; % ms
        mySinc_dur = floor(mySinc_dur/dt/Nsegs) * Nsegs*dt; % ms
        mySinc = sinc(linspace(-1, 1,  mySinc_dur/dt))' ;
        mySinc_area = (sum(mySinc) * dt);
        mySinc =  B1_area180 * mySinc / mySinc_area;
        B1pulse = mySinc;
        B1seg_len = round(mySinc_dur /dt/Nsegs);
       
    case 'sinc_mod'
        
        mySinc_dur = 2; % ms
        mySinc_dur = floor(mySinc_dur/dt/Nsegs) * Nsegs*dt; % ms
        mySinc = sinc(linspace(-1, 1,  mySinc_dur/dt))' ;
        mySinc_area = (sum(mySinc) * dt);
        mySinc =  B1_area180 * mySinc / mySinc_area;
        B1pulse = mySinc;
        B1seg_len = round(mySinc_dur /dt/Nsegs);
        
    case 'BIR'
        myBIR_dur = 2; % ms
        myBIRpulse = genBIRpulse(1 , myBIR_dur)';
        myBIRpulse_area = (sum(myBIRpulse) * dt);
        myBIRpulse = 3* B1_area180 *  myBIRpulse / myBIRpulse_area ;
        B1pulse = myBIRpulse;
        B1seg_len = round(myBIR_dur/dt/Nsegs);
        
    case 'sech'
        
        mysech_dur = 4 ;
        mysech = genSech180(4, mysech_dur)';
        mysech_area = sum(mysech * dt);
        mysech =  B1_area180 * mysech / mysech_area;
        B1pulse = mysech;
        B1seg_len = round( mysech_dur /dt/Nsegs);
        
    case 'hard'
        
        hard180_dur = 1.44; % (ms) this makes it an even multiple of 4 us and 9 segments
        hard180_dur = 0.16 * 9 ; % (ms) in the paper
        
        hard180_dur = round(hard180_dur/Nsegs/dt) * Nsegs*dt;
        hard180 = ones(hard180_dur/dt);
        
        hard180_area = sum(hard180)*dt;
        hard180 = B1_area180* hard180 /hard180_area;
        B1pulse = hard180;        
        B1seg_len = round(hard180_dur /dt/Nsegs);

end


% refocuser pulse is 180 hard pulse for 0.8 ms
B1ref_dur = 0.8;   % ms
B1ref_dur = 1.0;   % ms

B1ref_len = B1ref_dur/dt;  % points
B1ref_max = (1/B1ref_dur)* 0.1175e-4; % Tesla ... amplitude of 1 ms hard 180
B1ref = B1ref_max * ones(B1ref_len,1) ; % * exp(i*pi/2);
B1ref_area = sum(abs(B1ref));


%B1ref(:) = 0;


% removing the double refusing pulses
%B1ref(:) = 0;

%{
% 7 lobe sinc 180 for 5 ms
B1sinc7 = sinc(linspace(-8,8, 5/dt));
hw = hamming(5/dt)';
B1sinc7 = B1sinc7 .* hw;
B1sinc7_area = sum((B1sinc7));
B1sinc7 = B1sinc7 * B1ref_area / B1sinc7_area ;
B1sinc7_max = max((B1sinc7));
B1sinc7_area = sum((B1sinc7));
plot(B1sinc7)

sinc7_scale = B1sinc7_max/B1ref_max

% 1 lobe sinc 180 for 3.2 ms
B1sinc1 = sinc(linspace(-4,4, 3.2/dt));
hw = hamming(3.2/dt)';
B1sinc1 = B1sinc1 .* hw;
B1sinc1_area = sum((B1sinc1));
B1sinc1 = B1sinc1 * B1ref_area / B1sinc1_area ;
B1sinc1_max = max((B1sinc1));
B1sinc1_area = sum((B1sinc1));

plot(B1sinc1)
sinc1_scale = B1sinc1_max/B1ref_max

%}

str = sprintf('Durations: B1seg= %0.2f  ms, B1ref=  %0.2f ms , Gtrap=%0.2f ms\n', ...
    B1seg_len*dt , B1ref_len*dt, trap_len*dt);
%%
Gz = [];
B1 = [];

B1seg = B1pulse;

switch refocus_scheme
    case 'MLEV'
        % The MLEV-16 sequence
        mlev_phases =  repmat( [1 1 -1 -1  -1  1 1 -1 -1 -1  1 1  1  -1 -1  1], [1,4]);
        mlev_phases =  repmat( [1 1 -1 -1  -1  1 1 -1 -1 -1  1 1  1  -1 -1  -1], [1,4]);  % trying something different?
        
    case 'DR180'
        % My scheme for refocusing the off-res.
        mlev_phases = repmat( exp(-i*pi/2)*[1 -1 -1 1], [1,16]);
end

for n=1:Nsegs-1
    B1seg = B1pulse(1 + B1seg_len*(n-1) : n*B1seg_len);
    Gz = [Gz;
        zeros(size(B1seg)) ;
        trap;
        zeros(size(B1ref))
        -trap;
        trap;
        zeros(size(B1ref))
        -trap;
        ];
    
    if isControl
        Gz = abs(Gz);
        %Gz(:)  = 0;
    end
    
    
    B1 = [B1 ;
        B1seg;
        zeros(size(trap));
        B1ref * mlev_phases(2*n-1);
        zeros(size(trap));
        zeros(size(trap));
        B1ref * mlev_phases(2*n) ;
        zeros(size(trap));
        ];
end

Gz = [
    Gz;
    zeros(size(B1seg))
    ];

% control case:
if isControl
    Gz = abs(Gz);
end

B1seg = B1pulse(1 + B1seg_len*(Nsegs-1) : Nsegs*B1seg_len);

B1 = [B1;
    B1seg;
    ]*B1fudge;

if exportPulses
% write pulses to put on the scanner
    genScannerPulses(B1, Gz, dt);
end

Bx = real(B1);
By = imag(B1);

NSTEPS = length(B1);


% total duration of the simulation interval (in ms)
duration = NSTEPS*dt;  % ms.
vel_range =[-100:100]*1e-3;
vel_range =[-10:0.01:10]*1e-3;

for vel = vel_range  % cm / msec
    
    t = linspace(0,duration, NSTEPS)'; % mseconds.
    zpos = pos0 + vel*t + 0.5*accel*(t.^2);
    Bz = zpos.*Gz ;
    
    Bz = Bz + off_resonance;
        
    beff = [Bx By  Bz];
    Mi = [0 0 1]';
    M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
    if showEvolution
        plot(M(:,3))
        drawnow
    end
    
    Mzfinal=[Mzfinal; M(end,3)];
end


% now do the control case

for vel = vel_range  % cm / msec
    
    Bz = zpos.*abs(Gz) ;

    Bz = Bz + off_resonance;
        
    beff = [Bx By  Bz];
    Mi = [0 0 1]';
    M = blochsim(Mi, beff, T1, T2, dt, NSTEPS);
    if showEvolution
        plot(M(:,3))
        drawnow
    end
    
    Mzfinal_ns=[Mzfinal_ns; M(end,3)];
end


if batch_mode==0
%%    figure(3)
    t=[0:length(Gz)-1]*dt;
    
    subplot(311)
    area(t, Gz *1e4);
    grid on
    xlabel('time (ms)')
    ylabel('G_z (G/cm)')
    title('VSI Pulse')
    
    subplot(312)
    area(t, abs(B1)*1e7);
    grid on
    xlabel('time (ms)')
    ylabel('B_1 amplitude (mG)')
    
    subplot(313)
    area(t, angle(B1));
    grid on
    xlabel('time (ms)')
    ylabel('B_1 phase (rad)')
    
    figure
    hold on
    plot(vel_range*1e3, Mzfinal)
    plot(vel_range*1e3, Mzfinal_ns)
    fatlines
    grid on
    xlabel('Velocity (cm/s)')
    ylabel('M_z')
    title('Velocity Profile of VSI pulse')
    legend('Selective Case', 'Non-selective Case', 'Location', 'SouthEast')
    
    %text(-90 , -0.6, {refocus_scheme ;base_pulse });
    hold off
    
%%
end

%%  Calculating labeling efficency over a range of values
%{
figure(8)
efficiency = (Mzfinal - Mzfinal_ns)/2;
plot(vel_range, efficiency)
inds = find((vel_range < 0.08)  & (vel_range > 0.005));
hold on
plot(vel_range(inds), efficiency(inds))
mean(efficiency(inds))
xlabel('Velocity (cm/s)')
ylabel('Efficiency')
title('Velocity Profile of VSI pulse')
legend('Selective Case', 'Non-selective Case', 'Location', 'SouthEast')
hold off

%}


