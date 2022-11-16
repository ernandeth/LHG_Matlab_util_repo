% optimization of parameters for concentric shield around a figure-eight
% coil

close all
global FOV Nvox
global hst
global E1 Ex1 Ey1 Ez1
FOV=0.24;
Nvox=33;

% Calculate E field from TMS coil alone:
current =1e3/1e-4; % this is actually dI/dt in Amps/sec. realistic for TMS
[E1 Ex1 Ey1 Ez1] =  make_fig8(current, 0.03, 0, 0.085, 0);
save defaultFig8.mat E1 Ex1 Ey1 Ez1

% load up parameters for default coil from a file.  We do this only once
load defaultFig8.mat


hst.penetration=[];
hst.sharp=[];
hst.parms = [];
hst.cost=[];
hst.Ecylmax = [];
hst.W=1;


% init. paramaters for optiimzation 
sfactor = 0;
xgap = 0.02;
theta= 0 ;
R = 0.05;
Z = 0.25;
consts = [];

parms_0 = [
    sfactor   R    Z ] ;

LB = [  -3,  0.035,  0.09];
UB = [  3,   0.12,   0.25];

% Use these if you want to place the shield below the main coil"
%LB = [  -3,  0.035,  0.076];
%UB = [  3,   0.1,   0.084];

% Try a bunch of penetration weights
count=1; 
for W=[0 2 3 4:2:12]
    hst.penetration=[];
    hst.sharp=[];
    hst.parms = [];
    hst.cost=[];
    hst.Ecylmax = [];
    % hst.W=1;
    hst.W=W;
    
    [parms fval exitflag output] = fmincon(@shieldedTMS_lsq2, parms_0,[],[],[],[], LB, UB);
    
    allhst{count}=hst;
    count = count+1
    
    save bestParms_temp.mat allhst
    figure(45)
    subplot(221)
    plot(hst.sharp); title('Sharpness')
    subplot(222)
    plot(hst.penetration); title('Penetration')
    subplot(223)
    plot(hst.parms); title('Parameters'); legend('S.factor', 'Radius', 'Z-position'); 
    subplot(224)
    plot(hst.cost); title('Cost')

end


%%
