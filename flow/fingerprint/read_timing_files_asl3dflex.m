function [parms phys_parms] = read_timing_files_asl3dflex(myDir)
% function parms = read_timing_files(myDir)
%
% read acquisition parameter files for Velocity Selective ASL
% Fingerprinting time series experiment
% 
curDir = pwd;
cd(myDir)

parms.del1 =        1e-6*load('tadjusttbl.txt');
parms.del2 =        1e-6*load('prep1_pldtbl.txt');
parms.del3 =        1e-6*load('prep2_pldtbl.txt'); % delay between AS pulse and acqusition
parms.labelcontrol= load('prep1_lbltbl.txt');
parms.doArtSup =    load('prep2_lbltbl.txt');
parms.order =       load('order.txt');
parms.RO_time=      1e-6*load('RO_time.txt');    
parms.doblksat =    load('doblksattbl.txt'); 

parms.t_aq =        1e-6*load('RO_time.txt');%load('t_aqs.txt');
% hard code these:
fp = fopen('RO_type.txt','r')
parms.RO_type =     fscanf(fp,'%s',1);
fclose(fp);
fp = fopen('label_type.txt','r')
parms.label_type =     fscanf(fp,'%s',1);
fclose(fp);


% some nice defaults for simulation
    phys_parms.f = 0.01;
    phys_parms.Mtis0 = 1;
    phys_parms.cbva = 0.02;
    phys_parms.bat =  0.07;
    phys_parms.r1tis =  1/1.3;
    phys_parms.flip =  pi;
    phys_parms.r2tis = 1/0.090;

    total_dutaion = sum(parms.del1 + parms.del2 + parms.del3 +parms.RO_time)

cd(curDir)
return