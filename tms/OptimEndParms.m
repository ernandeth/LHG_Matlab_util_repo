% clear all
% close all
load bestParms_temp

% weight of penetration in cost function
Nsims = length(allhst)
for n=1:Nsims
    W(n) = allhst{n}.W;
    S(n)  = allhst{n}.sharp(end);
    P(n) = allhst{n}.penetration(end);
    F(n) = allhst{n}.parms(end,1);
    R(n) = allhst{n}.parms(end,2);
    Z(n) = allhst{n}.parms(end,3);
end
S0  = allhst{1}.sharp(1);
P0 = allhst{1}.penetration(1);
F0 = allhst{1}.parms(1,1);
R0 = allhst{1}.parms(1,2);


figure

plot(W,100*(S-S0)/S0);
hold on
plot(W, 100*(P-P0)/P0, '--');
hold off

legend('Sharpness','Penetration')
xlabel('Relative Weight of Penetration')
ylabel('% Change')
title('Effect of penetration weighting')
%fatlines; 
axis([0 10 -100 100])
dofontsize(16);
grid off

startparms = [S0 P0]
optparms = [W; S; P; F; R; Z]'
optparms = [1:Nsims;  W; 100*(S-S0)/S0; 100*(P-P0)/P0; 100*F; 100*R; 100*Z]'



