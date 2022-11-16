% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time


%%
% this chunk is for calculating inversion efficiency as a function of
% location of the inversion tag

doBatch=1;
sigs=[];
for tag_loc = -4:.5:8
% for RF_theta = -pi:.1:pi
%     tag_loc = 0;
    %RF_theta=0;
    flip_ang = 25;
    isTag = 1;
    vel = 30;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.4;  % this is the gradient imbalance used to gain some phase
    zpos = 0;
    
    pCASL01


    tagSignal = M(end,3);
    sigs = [sigs; tag_loc tagSignal Gz_phaser];
end

plot(sigs(:,1) , sigs(:,2))
plot(sigs(:,1) , sigs(:,3) ,'r')
plot(sigs(:,3) , sigs(:,2))

pause

isTag = 0;
pCASL01
conSignal = M(end,3);

alpha  = 0.5*(conSignal - tagSignal)/conSignal

%% this chunk is for stationary spins in different slices.
% requires commenting out the dz part of pCASL01

sigs=[];
doBatch=1;
for zpos = -5:0.5:5
    
    flip_ang = 25;
    
    vel = 0;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    tag_loc = 3.5 ; % cm
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.36;  % this is the gradient imbalance used to gain some phase 
    
    isTag = 1;
    pCASL01
    tagSignal = M(end,3);
    zpos
    
    
    % reset all the variables before running again.
    flip_ang = 25;
    vel = 0;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    tag_loc = 3.5 ; % cm
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.36;  % this is the gradient imbalance used to gain some phase
    
    isTag = 0; 
    %pCASL01
    conSignal = M(end,3);
    
    
    sigs = [sigs ;
        tagSignal conSignal phi];
  
end

figure
plot(-5:0.5:5, sigs(:,1)-sigs(:,2))
plot(-5:0.5:5, sigs(:,2))

%% this chunk is to simulate the effect of different degrees of gradient
% imbalance on the stationary spins.

sigs=[];
doBatch = 1;
for Gz_phaser = -2:.2:2

    flip_ang = 25;
    isTag = 1;
    vel = 30;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    tag_loc = 0 ; % cm
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    
    pCASL01
    tagSignal = M(end,3);
    
%    isTag = 0;
    
%    pCASL01
    conSignal = M(end,3);
    
    
    sigs = [sigs ;
        tagSignal conSignal Gz_phaser];

end

%plot(-1:.02:0.1,sigs(:,1))
plot(sigs(:,3),sigs(:,1))

%%
% this chunk is for calculating inversion efficiency as a function of
% velocity of the spins

doBatch=1;
sigs=[];
for velocity = [0:1:10 15 25 30]
    RF_theta2 = pi/4;

    
    tag_loc = 0;
    flip_ang = 25;
    isTag = 1;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.4;  % this is the gradient imbalance used to gain some phase
    zpos = 0;
    
    pCASL01
    tagSignal = M(end,3);

    tag_loc = 0;
    flip_ang = 25;
    isTag = 0;
    RF_spacing = 0.920 ; % time between RF pulses in ms.
    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    Gz_phaser = -0.4;  % this is the gradient imbalance used to gain some phase
    zpos = 0;
    
    pCASL01
    conSignal = M(end,3);
    
    sigs = [sigs; velocity tagSignal conSignal aq_phase];
    
end

figure
plot(sigs(:,1) , ( sigs(:,3)-sigs(:,2) )/2, 'k' )
grid on

%% Make some figures from saved data
load flipAngVelocity

figure
hold on
plot(sigs(:,1) , ( sigs10(:,3)-sigs10(:,2) )/2, 'k' )
plot(sigs(:,1) , ( sigs15(:,3)-sigs15(:,2) )/2, '.k' )
plot(sigs(:,1) , ( sigs20(:,3)-sigs20(:,2) )/2, '*k' )
plot(sigs(:,1) , ( sigs25(:,3)-sigs25(:,2) )/2, '--k' )
plot(sigs(:,1) , ( sigs35(:,3)-sigs35(:,2) )/2, '.-k' )

legend('\theta = 10', '\theta = 15', '\theta = 20', '\theta = 25', '\theta = 35')
xlabel('Velocity (cm/s)')
ylabel('Inversion efficiency ');
grid on
dofontsize(16)

figure
hold on
plot(sigs(:,1) , ( sigs_mpi4(:,3)-sigs_mpi4(:,2) )/2, '--k' )
plot(sigs(:,1) , ( sigs_0(:,3)-sigs_0(:,2) )/2, 'k' )
plot(sigs(:,1) , ( sigs_pi4(:,3)-sigs_pi4(:,2) )/2, '.-k' )
xlabel('Velocity (cm/s)')
ylabel('Inversion efficiency ');
legend('\phi_{extra} = -\pi/4' , '\phi_{extra} = 0' , '\phi_{extra} = \pi/4' )
grid on
dofontsize(16)

plot(sigs(:,1) , sigs(:,3))

%%