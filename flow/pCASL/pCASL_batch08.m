% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time


doBatch=1;
sigs=[];
phis = [];

for phase_correction = 0:0.1:2*pi
    tag_loc = 4;   
    flip_ang = 22.5;
    isTag = 1;
    vel = 50;
    slomo =0;
    eta = 0.35;   %
    t_ramp = 0.02;  % ramp time in ms.
    Gz_err = 0.05;  % artefactual gradient.  accounting for susceptibility
    

    % eta = 0.5 means balanced slice-select and refocuser : 
    % the refocuser area is half the ss area
    slomo = 0;
   
    pCASL08

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