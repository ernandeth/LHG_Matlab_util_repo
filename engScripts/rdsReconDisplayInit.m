%  This is the Initialization of the Recon, display, subtraction loop
%  in the Real Time ASL sequence.
%
%  first thing we do is figure out the dimensions from a sample P file
%  and calculate the kspace trajectory.  These things get used over and over and over and over
%

close all
global RT_MODE


RT_MODE=1;

DEBUG = 1;
DOCORRELATION=0;
DOREGRESSION=1;
SAVEMOVIE = 0;
SAVEMASK=0;
SMOOTHFLAG=1;
DEMEANFLAG=0;
DO3D = 0;
DORECON = 0;
NROWS=2;
ASLORDER=1;

InitialROIset=0;
myAxis=0;

figure(1);
set(gcf,'Position',[50 500 600 300], 'Toolbar', 'none', 'Menubar', 'none','Name', 'Last raw. image')
colormap(gray(256)); axis xy

figure(2);
set(gcf,'Position',[650 500 600 300],  'Toolbar','none', 'Menubar', 'none', 'Name', 'Last sub. image')
colormap(gray(256));axis xy

figure(3);
set(gcf,'Position',[50 100 600 300],  'Toolbar','none', 'Menubar', 'none', 'Name','Cummulative Map')
mymap = rt_make_colormap;
colormap(mymap);axis xy

figure(4);
set(gcf,'Position',[650 100 600 300],  'Toolbar','none', 'Menubar', 'none','Name','Time courses')


if exist('ReconParms.mat')
    fprintf('\nReading ReconParms.mat ... ');
    load ReconParms
    
elseif exist('Psample')
    % read header info if available
    fprintf('\nReading Psample ... ');
    [args,scaninfo,kinfo] = rec_setup1('Psample');
    Nslices = scaninfo.nslices;
    Nframes = scaninfo.nphases;
    disdaq = scaninfo.dda;
    xres = scaninfo.opxres;
    
else
    fprintf('No Parms file found ... use defaults ')
    
    scaninfo.nslices = 20;
    scaninfo.nphases = 200;
    scaninfo.dda = 2;
    scaninfo.opxres = 64;
end




save ReconParms_new.mat args scaninfo kinfo

% setup variables for the display loop.
xres = scaninfo.opxres;
Nframes = scaninfo.nphases;
Nslices = scaninfo.nslices;

Nrows = floor(sqrt(Nslices));
Nrows = 3;
Ncols = ceil(Nslices/Nrows);

ndat = scaninfo.ndat ;
Npix = xres*xres*Nslices;
var_con = zeros(1,Npix);

allRawData = zeros(Nframes, Nslices*scaninfo.ndat*2);
volumes = single(zeros(Nframes, Nslices*xres*xres));
allASL = single(zeros(Nframes , xres*xres*Nslices));

xc = 1;
yc = 1;

if DEBUG
    if 0
        [raw h] = read_img('nback.nii');
        Nframes = h.tdim;
        xres = h.xdim;
        Nslices = h.zdim;
        Npix = Nslices * xres *xres;
    end
    % alternatively we can do this:
    load RealTimeASL
    
    % override flags that we read in from file:
    %     DEBUG = 1;
    %     SAVEMOVIE=1;
    %     DOCORRELATION=0;
    %     SAVEMASK=0;
    %     SMOOTHFLAG=0;
    %     DEMEANFLAG=0;
    %     InitialROIset=0;
    %
    
    raw = volumes;
    clear allASL
    
end
% allocate space:
slicenum = 1;
vnum = 1;

begin = 1;
finish = xres*xres;

residual = single(zeros(Nframes , xres*xres*Nslices));
clean = single(zeros(Nframes , xres*xres*Nslices));
cumasl = single(zeros(1 , xres*xres*Nslices));
var_est = single(zeros(1 , xres*xres*Nslices));
var_con = single(zeros(1 , xres*xres*Nslices));
corrMap = single(zeros(1 , xres*xres*Nslices));
slmask = single(zeros(1 , xres*xres*Nslices));
tcrunch = zeros(Nframes,1);

mean_asl_course = zeros(Nframes,1);
wmjunk=0;
wm_inds = 1;
ref_inds = 1;
% use a variance mask to exclude high variance pixels from analysis

if exist('varMask.mat')
    load varMask
    SAVEMASK=0;
else
    varMask=ones(1,xres*xres*Nslices);
end

%Nframes = 200
if DOCORRELATION
    load('ReferenceWave.mat');
    
    ssxx=0;
    ssyy=0;
    ssxy=0;
    ssx2=0;
    ssy2=0;
    
    corr_map = zeros(size(varMask));
end
if DOREGRESSION
    load('ReferenceWave.mat');
end
% in debugging mode, loop through the recon and display stuff
% otherwise, let the rdsClient initiate the recon/display

if SAVEMOVIE
    movavi=avifile('RTmovie.avi','FPS',5,'COMPRESSION','MSVC','QUALITY',100);
end



order = 1:Nslices;

%datacursormode on
if DO3D
    % allocate space for a 3D FFT buffer
    tmp3d = single(zeros(xres,xres, Nslices ));
    tmp3drec = single(zeros(xres,xres, Nslices ));
    
    
    for n=1:Nslices/2
        order(2*n-1) = Nslices/2 -n+1;
        order(2*n) = Nslices/2 +n;
    end
    
    
end

if DEBUG
    profile on
    
    
    for c=1:Nslices*Nframes
        rdsReconDisplayLoop01
    end
    
    profile viewer
end

