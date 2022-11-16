% optimization of parameters for concentric shield around a figure-eight
% coil

close all
global FOV Nvox
global hst
global E1 Ex1 Ey1 Ez1
FOV=0.24;
Nvox=33;

% load up parameters for default coil from a file.  We do this only once
load defaultFig8.mat



%% init. paramaters for optiimzation 
sfactor = -0.7;
R = 0.040;
Z = 0.095;

hst.penetration=[];
hst.sharp=[];
hst.parms = [];
hst.cost=[];
hst.Ecylmax = [];
hst.W=1;

for sfactor= linspace(-3,3,20)
    parms_0 = [ sfactor   R    Z ] ;
    [E cost] = shieldedTMS_lsq2( parms_0,[]);
    

    count = count+1;
end
    allhst{1}=hst;
figure(33)
subplot(1,3,1)

plot(linspace(-3,3,20) , hst.penetration)
hold on
plot(linspace(-3,3,20) , 10*hst.sharp, '--')
axis tight
xlabel('Shield Current Factor')
%% init. paramaters for optiimzation 
sfactor = -0.7;
R = 0.040;
Z = 0.095;
hst.penetration=[];
hst.sharp=[];
hst.parms = [];
hst.cost=[];
hst.Ecylmax = [];
hst.W=1;
 
for R= linspace(0.035,0.12,20)
    parms_0 = [ sfactor   R    Z ] ;
    [E cost] = shieldedTMS_lsq2( parms_0,[]);

end

subplot(1,3,2)

plot(linspace(3.5,12,20) , hst.penetration)
hold on
plot(linspace(3.5,12,20) , 10*hst.sharp, '--')
legend( 'Penetration','Sharpness','Location','NorthWest')
legend boxoff
axis tight
xlabel('Shield Radius (cm)')
drawnow
%% init. paramaters for optiimzation 
sfactor = -0.7;
R = 0.040;
Z = 0.095;
hst.penetration=[];
hst.sharp=[];
hst.parms = [];
hst.cost=[];
hst.Ecylmax = [];
hst.W=1;
 
for Z= linspace(0.09, 0.25, 20)
    parms_0 = [ sfactor   R    Z ] ;
    [E cost] = shieldedTMS_lsq2( parms_0,[]);
end
allhst{3}=hst;
subplot(1,3,3)

plot(linspace(9, 25, 20) , hst.penetration)
hold on
plot(linspace(9, 25, 20) , 10*hst.sharp, '--')
xlabel('Shield Distance (cm)')
axis tight
%%
