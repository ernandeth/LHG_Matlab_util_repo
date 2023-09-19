function [parms phys_parms] = read_timing_files(myDir)
% function parms = read_timing_files(myDir)
%
% read acquisition parameter files for Velocity Selective ASL
% Fingerprinting time series experiment
% 
curDir = pwd;
cd(myDir)

parms.t_tags =      load('t_tags.txt');
parms.del1 =        load('t_adjusts.txt');
parms.del2 =        load('t_delays.txt');
parms.del3 =        load('AS_delays.txt'); % delay between AS pulse and acqusition
parms.labelcontrol= load('isVelocitySelective.txt');
parms.doArtSup =    load('doArtSuppression.txt');
parms.order =       load('order.txt');
parms.RO_time=      load('RO_time.txt');    


parms.t_aq =        load('RO_time.txt');%load('t_aqs.txt');
% hard code these:
parms.RO_type =     char(load('RO_type.txt'));
parms.label_type =  'BIR8inv'; % load('LabelType.txt');

% some nice defaults for simulation
    phys_parms.f = 0.01;
    phys_parms.Mtis0 = 1;
    phys_parms.cbva = 0.02;
    phys_parms.bat =  0.07;
    phys_parms.r1tis =  1/1.3;
    phys_parms.flip =  pi;
    phys_parms.r2tis = 1/0.090;

    total_dutaion = sum(parms.del1 + parms.del2 + parms.del3 +parms.RO_time+parms.t_tags)

cd(curDir)
return