% Batch file for pseudoCASL simulations.  COntains a bunch of cases
% execute a chunk of code at a time

% June 20, 2008
% this batch checks velocity dependence using different RF_pulsewidt
% params.

doBatch=1;
sigs=[];
phis = [];
vel = 20;


for Gz_ss = 0:0.5:3
    
    tag_loc = 0;   
    flip_ang = 22.5;
    isTag = 1;
    vel = 10;
    slomo =0;
    eta = 0.25;   % works for all locs if eta = 0.25!!
    t_ramp = 0.02;  % ramp time in ms.
    Gz_err = 0;  % artefactual gradient.  accounting for susceptibility
    RF_spacing = 1.5 ;  % default = 0.8.  check below.
    vel_phase = 0; 
    
    pCASL05

    tagSignal = (M(end,3));
    sigs = [sigs; tag_loc tagSignal];

end

plot(  sigs(:,2))
axis([0 40 -1 1])
grid on

return

% 
% parms that work:
%     flip_ang = 22.5;
%     isTag = 1;
%     vel = 40;
%     slomo =0;
%     eta = 0.25;   % works for all locs if eta = 0.25
