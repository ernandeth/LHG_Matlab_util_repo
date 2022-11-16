% 
% ObjectiveFunction = @gmin_mix1;
% nvars = 60;    % Number of variables
% random_l(1:20) = linspace(1,1,20);
% random_u(1:20) = linspace(1,1,20)*3;
% random_l(21:60) = linspace(0,0,40);
% random_u(21:60) = linspace(0,0,40)*1.5;
% 
% LB = random_l;   % Lower bound
% UB = random_u;  % Upper bound
% rng(1,'twister'); % for reproducibility
% [x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectiveFunction = @gmin_mix1;
nvars = 60;    % Number of variables
random_l(1:20) = linspace(1,1,20)*2.4;
random_u(1:20) = linspace(1,1,20)*5;  %%% TR is from 2.4 to 5sec
random_l(21:60) = linspace(0,0,40);
random_u(21:60) = linspace(0,0,40)*0.85;  %%% PID and lable_dur is from 0 to 0.85 sec

LB = random_l;   % Lower bound
UB = random_u;  % Upper bound
rng(1,'twister'); % for reproducibility
[x,fval] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB);