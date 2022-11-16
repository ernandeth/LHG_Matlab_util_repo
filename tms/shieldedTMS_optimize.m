close all

optvar = optimset('lsqnonlin');
%optvar = optimset('fminimax');
optvar.Display='on';
%optvar.TolFun = 1e-15;
%optvar.TolX = 1e-10;
%optvar.MaxIter = 100;
optvar.Diagnostics = 'on';
optvar.DiffMinChange = 1.0e-6;
%optvar.Algorithm='trust-region-reflective';

% init. paramaters for optiimzation 
sfactorBottom =0.5;
xgap = 0.02;
theta= pi/8 ;
Rbottom = 0.07;
%Rtop = 0.03;    
Zbottom = 0.15;
%top_Z = 0.15;
consts = [];

parms_0 = [
    sfactorBottom    xgap    theta    Rbottom    Zbottom
    ] ;

% % case 1: top shield only
% parms_0 = [ 0.5   0.01  pi/4  0.0300    0.20];
% 
% LB = [     -2, -0.15,  -pi/4,    0.03,    0.13 ];
% UB = [     2,   0.15,  pi/4,     0.10,    0.20   ];
% aperture = pi/8;

% case 2:  bottom shield only
%parms_0= [ -0.5   0.05   -pi/4    0.1000    0.1100 ];
%LB = [    -2, -0.05,  -pi/4,    0.07,    -0.13    ];
%UB = [    2,   0.05,  pi/4,     0.12,    0.12	];

% case 3:   shield between fig 8 and head sphere: - these give good results
% when optimizing on the E field!
% parms_0= [ -0.5   0.05   -pi/4    0.07000    0.1100 ];
% LB = [    -2, -0.05,  -pi/4,    0.02,    0.08    ];
% UB = [    2,   0.05,  pi/4,     0.06,    0.125	];

% case 3:   shield between fig 8 and head sphere:
parms_0= [ -0.5   0.05   -pi/4    0.04000    0.1300 ];
LB = [    -2, -0.05,  -pi/4,    0.02,    0.125    ];
UB = [    2,   0.05,  pi/4,     0.06,    0.25	];

%[parms resnorm residual exitflag output] = lsqnonlin(@shieldedTMS_lsq, parms_0, LB, UB, optvar); 
% pen=22, sharp=.12 boundaries were violated

% [parms fval] = fgoalattain(@shieldedTMS_lsq, parms_0,[0 0],[1 1],[],[],[],[], LB, UB);  
% pen=22, sharp=0.06
% fgoalattain never changed the parms!!

[parms fval] = fminimax(@shieldedTMS_lsq, parms_0,[],[],[],[], LB, UB); 
% pen=22,  sharp=0.13

sfactorBottom = parms(1);
%sfactorTop = parms(2);
xgap = parms(2);
theta = parms(3);
%Rtop = parms(4);
Rbottom = parms(4);
Zbottom = parms(5);
%top_Z = parms(7);

save bestParms_smallBottom_div.mat 
parms_baseline = [ 0 0 0 0 1];
%%
