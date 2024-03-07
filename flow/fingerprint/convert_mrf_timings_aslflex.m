% convert parameter files to asl3dflex format
% units of of the timing fiels are in seconds
!mkdir scnr_timing

% read parms in current format
aqparms = read_timing_files('./')
Nframes = length(aqparms.del1);
% write parms in asl3dflex format

% prep pulse 1i type of pulse
p1_id = aqparms.label_type;

switch(p1_id)
    case 'BIR8inv'
        p1_id=6800;
    case 'BIR8'
        p1_id=6850;round(tmp*1e6)
    case 'FTVSI-sinc'
        p1_id=17268;
end    
wtextfile('scnr_timing/prep1_id.txt', p1_id, '%d')

% prep pulse 2 type of pulse
p2_id = 6850;  % hardcode as a BIR8 pulse
wtextfile('scnr_timing/prep2_id.txt', p2_id, '%d')


% Global saturation : The first one is set to zero, otherwise it's always on
doGS = ones(Nframes,1);
doGS(1) = 0;
wtextfile('scnr_timing/doblksattbl.txt',  doGS, '%d')
GS_duration = doGS*2.5e-3; % force the global sat pulse to be 2.5 ms long.

% delay before first label pulse, aka, t_adjust or del1
% time between global sat and begining of label pulse 1
tmp = aqparms.del1 - GS_duration ;
tmp = round(tmp*1e6);  % from seconds to microseconds
wtextfile('scnr_timing/tadjusttbl.txt', tmp, '%d')

%-------------
% do nothing, label or control?
tmp = aqparms.labelcontrol;
wtextfile('scnr_timing/prep1_lbltbl.txt', tmp, '%d')
conve

% prep1 pld , aka, t_delay or del2
% time between ibeginning of label pulse 1 and beginging of label pulse 2
tmp = aqparms.del2 - p1_id*4e-6;
tmp = round(tmp*1e6);  % from seconds to microseconds
wtextfile('scnr_timing/prep1_pldtbl.txt', tmp, '%d')


%---------
% prep pulse 2: do nothing, label or control?
tmp = aqparms.doArtSup;
wtextfile('scnr_timing/prep2_lbltbl.txt',  tmp, '%d')

% prep2 pld: aka, ArtSup_delay, aka, del3
% subtract duration of art sup pulse, add the fatsat delay here
% N.B do not include fat sat time into the read out!!
tmp = aqparms.del3 - p2_id* 4e-6 - 7.5e-3;
tmp = round(tmp*1e6);  % from seconds to microseconds
wtextfile('scnr_timing/prep2_pldtbl.txt',  tmp, '%d')
