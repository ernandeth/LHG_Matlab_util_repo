%  This is the Initialization of the Recon, display, subtraction loop
%  in the Real Time ASL sequence.
%  
%  first thing we do is figure out the dimensions from a sample P file
%  and calculate the kspace trajectory.  These things get used over and over and over and over 
%   

DEBUG = 0


if exist('Psample')
    % read header info if available
    [args,scaninfo,kinfo] = rec_setup1('Psample');
    Nslices = scaninfo.nslices;
    Nframes = scaninfo.nphases;
    disdaq = scaninfo.dda;
    xres = scaninfo.opxres;
else
    fprintf('No Parms file found ... use defaults ')

    scaninfo.nslices = 5;
    scaninfo.nphases = 200;
    scaninfo.dda = 2;
    scaninfo.opxres = 64;
end

save reconParms.mat args scaninfo kinfo

fidsize=scaninfo.ndat;


if DEBUG
    [raw h] = read_img('nback.nii');
    Nframes = h.tdim;
    xres = h.xdim;
    Nslices = h.zdim;
    
end
Nframes = 200

slicenum = 1;
vnum = 1;

begin = 1;
finish = xres*xres;

volumes = zeros(Nframes, Nslices*xres*xres);
allASL = zeros(Nframes-2 , xres*xres*Nslices);
mean_asl_course = zeros(Nframes-2,1);


% in debugging mode, loop through the recon and display stuff
% otherwise, let the rdsClient initiate the recon/display 
if DEBUG
    for c=1:Nslices*Nframes
        rdsReconDisplayLoop01
    end
end
