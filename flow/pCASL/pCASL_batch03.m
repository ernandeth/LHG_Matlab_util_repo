% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time


%% this chunk is for calculating inversion efficiency
doBatch=1;
sigs=[];
for tag_loc = -4:.5:8
    % for RF_theta = -pi:.1:pi
    %     tag_loc = 0;
    %RF_theta=0;
    flip_ang = 22.5;
    isTag = 0;
    vel = 30;

    RF_spacing = 1.80 ; % time between RF pulses in ms.

    t_ramp = 0.02;  % ramp time in ms.
    Gz_ss = 0.6;  % slice select gradient in G/cm.
    t_phaser = 300 ; % us
    Gz_phaser = -0.36;  % this is the gradient imbalance used to gain some phase

    pCASL03


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
for zpos = -12:2:12

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
    pCASL01
    conSignal = M(end,3);


    sigs = [sigs ;
        tagSignal conSignal phi];

end

figure
plot(-12:2:12, sigs(:,1)-sigs(:,2))
plot(-12:2:12, sigs(:,1))

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