close all
global FOV Nvox
FOV=0.5
FOV = 0.24
Nvox=33
%optvar = optimset('lsqnonlin');
optvar = optimset('fmincon');
optvar.Display='on';
%optvar.TolFun = 1e-15;
%optvar.TolX = 1e-5;
%optvar.MaxIter = 100;
optvar.MaxFunEvals = '200*numberofvariables';
%optvar.Diagnostics = 'on';
optvar.DiffMinChange = 1.0e-4;
%optvar.PlotFcns=@optimplotresnorm;
%optvar.PlotFcns=@optimplotfval;
%optvar.Algorithm='trust-region-reflective';


% parms_0 = [
%    sfactorBottom      Rbottom    Zbottom ] ;

parms_0 = [ 0  0.150  0.13 ];

LB = [     -2,  0.02  0.13];
UB = [     2,   0.18  0.20];

% [parms resnorm residual exitflag output ] =	lsqnonlin(@shieldedTMS_lsq3, parms_0, LB, UB, optvar); 
% pen=22, sharp=.12 boundaries were violated

% [parms fval] = fgoalattain(@shieldedTMS_lsq, parms_0,[0 0],[1 1],[],[],[],[], LB, UB);  
% pen=22, sharp=0.06
% fgoalattain never changed the parms!!

%[parms fval] = fminimax(@shieldedTMS_lsq3, parms_0,[],[],[],[], LB, UB, [],optvar);
[parms fval] = fmincon(@shieldedTMS_lsq3, parms_0,[],[],[],[], LB, UB);
% pen=22,  sharp=0.13


save bestParms_square.mat 
parms_baseline = [ 0 0 ];
%%
