% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time


doBatch=1;
sigs=[];
phis = [];

for zpos = -5:0.1:-1
    tag_loc = 3;
    flip_ang = 22.5;
    isTag = 1;
    vel = 0;
    slomo =0;
    eta = 0.4;   % works for all locs if eta = 0.25!!
    t_ramp = 0.02;  % ramp time in ms.
   
    pCASL06

    tagSignal = (M(end,3));
    sigs = [sigs; tag_loc tagSignal];

end

plot(sigs(:,2))


return

% 
% parms that work:
%     flip_ang = 22.5;
%     isTag = 1;
%     vel = 40;
%     slomo =0;
%     eta = 0.25;   % works for all locs if eta = 0.25