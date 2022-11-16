function [args,gkseries4D,imseries4D,im_psf] = spiral3drecon_220412(varargin)
%function [args,gkseries4D,imseries4D,im_psf] = spiral3drecon(varargin)
%
% spiral3drecon
%
% Luis Hernandez-Garcia, David Frey @UM 2020
%
% 1- This script reads FIDs from the GE scanner acquired with a
% noncartesian 3d trajectory
% 2 - it reads the gradient and kspace trajectories and rotates them
% appropriately for each shot - USES rotation matrices sorted in text files.
% 3 - reconstructs the image by Gridding and by Jeff's NUFFT
%
% Variable Inputs (varargin):
%   myPfile: name of Pfile; can set to 'auto' to to use first Pfile in
%       directory (default is 'auto')
%   saveArgs: option to save recon configuration to a .mat file and load
%       later (default is 1). If done so, recon can be re-ran with same
%       arguments by passing 'lastargs' as first parameter
%   doNUFFT: option to use NUFFT (default is 1)
%   doDensComp: option to incorporate sampling density compensation
%       (default is 'auto': set by doNUFFT)
%       dcfN: Number of neighbors to use in density calculation (default
%           is 100)
%       dcfR: Exponent on density calculation (default is 2)
%   zoomFactor: option to zoom in on reconned image by adjusting delta k
%       (default is 1)
%   scaleFactor: option to scale data; pass 'auto' to use normalized and scaled
%       by DAQ factor (default is 'auto')
%   Ndel: sample delay (default is -3)
%   SHOWPIX: debug level, shows more images along the way with higher
%       levels (default is 0)
%       levels:
%           1: display ov of gridded recon and kspace at end of each frame
%               recon loop
%           2: display figures showing raw data, trajectory & density, psf
%           3: changes parfor recon loop to regular loop (0 worker parfor)
%               and displays figures generated within the loop, plus
%               additional "slow" figures such as T2filter
%   MaxFrame: number of frames to recon; can set to 'auto' to assign value
%       based on SHOWPIX (all frames for SHOWPIX<3, 1 for SHOWPIX>=3); can
%       force use of all frames by setting to 'all'. (default is 'auto')
%   doRots: use rotating spiral platter in calculating k-space traj; can
%       set to 'auto' to assign value based on psd name (1 if psdname is 
%       vsasl3dflex04, 0 otherwise). (default is 'auto')
%   useKrot: if using rotating spiral, option to read rotation matrices
%       from krotations.txt rather than calculating new ones; can set to
%       'auto' to assign value based on doRots (=doRots). (default is
%       'auto')
%   doYrot: if using a rotating spiral with Y-axis (psi) rotations (default
%       is 0)
%   rotAngle: slice rotation angle for rotating spiral (default is golden
%       angle = pi*(3-sqrt(5)) )
%   doPhaseCycle: option to include CGMP phase cycling (default is 0)
%   doPhaseDetrend: option to perform phase detrending (default is 0), if
%       1, can also pass the following parameters:
%       PDtype: method of phase detrending (lsqNav, lsqAll, or HPF)
%           (default is lsqNav)
%       PDspec: specification for PD method (max order if lsq, cutoff if
%           HPF) (default is 0 -- Must change to >0 for HPF to work)
%       PDunwrap: option to unwrap phase for PD calculations (default is 0)
%   doT2filter: option to perform T2 filtering (default is 0)
%   T2type: specify how T2 values act on data
%       1: Each echo train (frame/coil/leaf) has an individual T2 fit
%       2: Each frame/coil has an individual T2 fit
%       3: One T2 fit for the entire dataset
%
%   To add a new input variable:
%       1. Append "defaults" structure in Section A to include the new
%           variable and default value
%       2. Append the line "<input parm name> = args.<input parm name>"
%           to the end of section B
%       3. (Optional) Add a validator function in spiral3drecon_inputparser
%           by defining a logical function handle that returns true if the
%           input is of the right format, and false otherwise. If this step
%           is skipped, the validator will always return true by default.
%
% Outputs:
%   args: structure containing configuration of all input arguments
%   gkseries4D: 4D dataset containing uncompressed/unscaled k-space data
%   imseries4D: 4D dataset containing uncompressed/unscaled image data
%   im_psf: 3D dataset containing point-spread function

%% Section A: Setting default constants & recon configurations:
defaults = struct(...
    'myPfile', 'auto', ... % default myPfile is automatically set from directory
    'saveArgs', 1, ...
    'doNUFFT', 1, ...
    'doDensComp', 'auto', ... % default doDensComp is automatically set from doNUFFT
    'dcfN', 100, ...
    'dcfR', 2.5, ...
    'zoomFactor', 1, ...
    'scaleFactor', 'auto', ... % default scaleFactor is automatically set from DAQfactor
    'Ndel', -3, ...
    'SHOWPIX', 0, ...
    'MaxFrame', 'auto', ... % default MaxFrame is automatically set from SHOWPIX
    'doRots', 'auto', ... % default doRots is automatically set from psd name
    'useKrot', 'auto', ... % default useKrot is automatically set from doRots
    'doYrot', 0, ...
    'rotAngle', pi*(3-sqrt(5)), ...
    'doPhaseCycle', 0, ...
    'doPhaseDetrend', 0, ...
    'PDtype', 'lsqNav', ...
    'PDspec', 0, ...
    'PDunwrap', 0, ...
    'doT2filter', 0, ...
    'T2type', 3 ...
    );

% Defining important constants:
dt = 4e-6; %seconds
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

DAQfac = 2^15 - 1;
Nramp = 20; % Number of ramp points to throw out

%% Section B: Parsing through variable inputs:
if numel(varargin)<1
% Skip input parser and set to defaults:
    args = defaults;
elseif strcmpi(varargin{1},'lastargs')
% Load previous configuration and use as defaults in input parser:
    load lastreconargs.mat args
    varargin(1) = [];
    args = spiral3drecon_inputparser(varargin,args);
else
% Parse inputs:
    args = spiral3drecon_inputparser(varargin,defaults);
end

% Set auto parameters:
if strcmpi(args.myPfile,'auto')
    pfiles = dir('P*.7');
    args.myPfile = pfiles(1).name; % auto assign myPfile based on dir:
end

    % read Pfile now (need pfile name to do so, parms need scaninfo to be set):
    fprintf("\nLoading Pfile: %s\n", args.myPfile)
    [raw,scaninfo] = read_raw_3d(args.myPfile,0);

if strcmpi(args.doDensComp,'auto') && (args.doNUFFT)
    args.doDensComp = 1; % auto assign doDensComp based on doNUFFT
elseif strcmpi(args.doDensComp,'auto') && ~(args.doNUFFT)
    args.doDensComp = 0; % auto assign doDensComp based on doNUFFT
end

if strcmpi(args.MaxFrame,'auto') && (args.SHOWPIX>=3)
    args.MaxFrame = 1; % auto assign MaxFrame based on SHOWPIX:
elseif strcmpi(args.MaxFrame,'auto') || strcmpi(args.MaxFrame,'all')
    args.MaxFrame = scaninfo.nphases; % auto assign MaxFrame based on SHOWPIX:
end

if strcmpi(args.doRots,'auto') && strcmpi(scaninfo.psdname,'vsasl3dflex04')
    args.doRots = 1; % auto assign doRots based on psd name:
elseif strcmpi(args.doRots,'auto')
    args.doRots = 0; % auto assign doRots based on psd name:
end

if strcmpi(args.useKrot,'auto')
    args.useKrot = args.doRots; % auto assign useKrot based on doRots:
end
    
if args.saveArgs % save recon args for later
    save lastreconargs.mat args
end

% Assign variables from args:
myPfile = args.myPfile;
saveArgs = args.saveArgs;
doNUFFT = args.doNUFFT;
doDensComp = args.doDensComp;
dcfN = args.dcfN;
dcfR = args.dcfR;
zoomFactor = args.zoomFactor;
scaleFactor = args.scaleFactor;
Ndel = args.Ndel;
Ndel2 = Ndel;
SHOWPIX = args.SHOWPIX;
MaxFrame = args.MaxFrame;
doRots = args.doRots;
useKrot = args.useKrot;
doYrot = args.doYrot;
rotAngle = args.rotAngle;
doPhaseCycle = args.doPhaseCycle;
doPhaseDetrend = args.doPhaseDetrend;
PDtype = args.PDtype;
PDspec = args.PDspec;
PDunwrap = args.PDunwrap;
doT2filter = args.doT2filter;
T2type = args.T2type;

%% Section C: Parsing out pfile header information:
nslices = scaninfo.nslices;
nechoes = scaninfo.npr;
nframes = scaninfo.nphases;
ndat = scaninfo.ndat;
ncoils = scaninfo.ncoils;
nfids = nechoes*nframes*nslices;
fov = scaninfo.opfov;
MYDIM = scaninfo.opxres;
nleaves = scaninfo.npr;
slthick = scaninfo.slthick/10;
TE = scaninfo.te;
fov3 = [fov fov ~doRots*nslices*slthick + doRots*fov];
dim3 = [MYDIM MYDIM ~doRots*nslices + doRots*MYDIM];
%%%%%%%%%%%%%%%%%%%%%%%%
% KLUGE (zoom factor does not work properly with nufft so reset it to 1)
%   Need to fix later:
if doNUFFT, zoomFactor = 1; end
%%%%%%%%%%%%%%%%%%%%%%%%
dim3 = round(dim3/zoomFactor)+1; % implement zoom factor and add one to ensure center of kspace is included

% define output header:
h = define_nii_hdr();
    h.dim = [4 dim3 nframes 0 0 0];
    h.pixdim = [1 fov3./dim3 1 0 0 0];
    h.datatype = 4;
    h.bitsperpixel = 16;
ah = nii2avw_hdr(h);

if SHOWPIX>=2
    % Plot raw echo train:
    fraw = figure;
    plot(squeeze((abs(raw(1,:,:)))));
    xlim([1 size(raw,2)])
    title('Echo train: Frame 1, all coils');
    xlabel('Data points');
    drawnow
end

fprintf('Beginning 3D spiral recon for %s with following configurations:\n', myPfile)
fprintf('\t%s\n',spiral3drecon_structinfo(args))

if saveArgs
    % write info to history file:
    ![ ! -f history.txt ] && touch history.txt
    fhistory = fopen('history.txt','a');
    fprintf(fhistory,'Recon of %s, %s\nRecon configurations:\n',myPfile,datetime);
    fprintf(fhistory,'%s\t',spiral3drecon_structinfo(args));
    fprintf(fhistory,'\nScan info:\n');
    fprintf(fhistory,'%s\t',spiral3drecon_structinfo(scaninfo));
    fprintf(fhistory,'\n\n\n');
end

%% Section D: Processing trajectory:
% Determine Ncenter from spherical trajectory:
ktraj_sph = load('ktraj_sph.txt');
ktraj_sph = ktraj_sph(2:end-1,1);
% Find points at center of kspace:
centerpts = find(ktraj_sph(:)  == 0);
% Exclude points in outer quartiles representing beginning/end of spiral:
centerpts(centerpts<ndat/4) = []; centerpts(centerpts>ndat*3/4) = [];
ncenter = length(centerpts);
scaninfo.ncenter = ncenter;

% Configure k-space trajectory:
fprintf('\nLoading grads...\n');
g = load('grad.txt');
    glen = size(g,1);
    aqwindow = dt*glen;

% Calculate the k-space trajectory calculated from the stored 
% gradient waveformss, rather than assuming an ideal trajectory
k2 = gamma*cumsum(g)*dt;  % cm-1 units
k2 =  k2';

% Calculate the k-space boundaries of the readout:
Kxymax = MYDIM/fov/2;  % Kmax
Kzmax = doRots*Kxymax + ~doRots*1/(2*slthick);

% Resample the trajectory to match the data acquisition rate:
ks = zeros(3,ndat);
for n=1:3
    tmp = k2(n,:);
    ks(n,:) = interp1(linspace(0,1,length(tmp)) , tmp, linspace(0,1,ndat));
end

% Throw out ramp/center points if specified:
ks = ks(:,Nramp+1:end-Nramp) ;
ks(:,round(end/2)-round(ncenter/2):round(end/2)+round(ncenter/2)) = [];

% Configurations for processing trajectory:
ks_all = [];
ks_plate = (doRots*2-1)*ks;
if useKrot, krot = load('krotations.txt'); end
if ~doRots, kzsteps = load('kzsteps.txt'); end

% Translating platter to each shot/slice:
![ -f krotations_calc.txt ] && rm -f krotations_calc.txt
!touch krotations_calc.txt
krot_calc = fopen('krotations_calc.txt','w');
for leaf = 1:nleaves
    for slice = 1:nslices
        idx_plate = (leaf-1)*nslices + slice;
        
        % Calculate rotation matrix:
        xi = doRots*rotAngle*(slice-1 + (leaf-1)/nleaves);
        psi = xi*doYrot;
        phi = pi*(leaf-1)/nleaves;
        Ryxz = -[0 1 0; 1 0 0; 0 0 1]*genRM3(["y" "x" "z"], [psi xi phi]);
        
        % Print to krotations_calc.txt for debugging:
        fprintf(krot_calc, '%d\t', [slice-1 leaf-1]);
        fprintf(krot_calc, '%.6f\t', [xi phi]);
        fprintf(krot_calc, '%d\t', round(DAQfac*Ryxz(:)'));
        fprintf(krot_calc, '\n');
        
        if useKrot
            % Determine rotation angles from krotations.txt:
            xi_krot = krot(idx_plate,3);
            phi_krot = krot(idx_plate,4);
            
            RMSE_ang = sqrt(mean(([xi_krot phi_krot]-[xi phi]).^2, 'all'));
            if RMSE_ang >= 1e-2
                warning('Leaf %d, Slice %d: rotation angles from krotations.txt do not match calculated angles (RMSE: %.3f)\n', ...
                    leaf, slice, RMSE_ang);
            end
            
            % Determine rotation matrix from krotations.txt:
            Ryxz_krot = reshape(krot(idx_plate,end-8:end),3,3)'/DAQfac;
            
            RMSE_R = sqrt(mean((Ryxz-Ryxz_krot).^2,'all'));
            if RMSE_R >= 1e-2
                warning('Leaf %d, Slice %d: rotation matrix from krotations.txt does not match calculated matrix (RMSE: %.3f)\n\tProceeding with matrix from krotations.txt anyways...', ...
                    leaf, slice, RMSE_R);
            end
            
            Ryxz = Ryxz_krot;
            
        end
        
        if ~doRots
            % For SOS, shift slices along zaxis:
            ks_plate(3,:) = Kzmax-kzsteps(slice,3);
        end
        
        Ryxz = [ 1 0 0; 0 -1 0; 0 0 -1]*Ryxz;
        
        ks_all = [ks_all; (Ryxz*ks_plate)'];
    end
end

ks = ks_all;

if SHOWPIX>=2
    % Plot trajectories:
    ftraj = figure;
    subplot(2,3,[1 4])
        x = 1:size(g,1);
        spacing = max(abs(g(:)))*2;
        plot(x,g(:,1)+spacing,x,g(:,2),x,g(:,3)-spacing)
        [minx,maxx] = bounds(x); xlim([minx maxx]), xticks([])
        yticks(spacing*(-1:1)), yticklabels(["z" "y" "x"])
        title(sprintf("First platter\ngradient components"))
    subplot(2,3,[2 3 5 6])
        plot3(ks(:,1),ks(:,2),ks(:,3))
        xlabel("x"), ylabel("y")
        title("3D trajectory")
    sgtitle("K-space Trajectory")
    drawnow
end

fprintf('k-space trajectory ready\n');

%% Section E: Pre-recon Processing:
% Initializing image timeseries:
imseries4D = zeros([dim3 nframes]);
imseries = zeros([prod(dim3) nframes]);
gkseries4D = imseries4D;
gkseries = imseries;

% Creating cartesian grid for recon:
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kxymax,Kxymax,dim3(1)), ...
    linspace(-Kxymax,Kxymax,dim3(2)), ...
    linspace(-Kzmax,Kzmax,dim3(3)));

% Calculating dists without ramp/center pts:
ktraj_sph = ktraj_sph(Nramp+1:end-Nramp);
ktraj_sph(round(end/2)-round(ncenter/2):round(end/2)+round(ncenter/2)) = [];
dists = repmat(ktraj_sph,nslices*nechoes,1);

% Create k-space sampling density function:
DCs = squeeze(mean(raw(:,round(end/2)-round(ncenter/2):round(end/2)+round(ncenter/2),:),2));
if doDensComp
    kdens = spiralDCF(ks,dcfN,dcfR);
else
    kdens = ones(size(ks,1),1);
end

% Form kdens into a grid and show image if specified:
kdens_grid = spiralGrid(ks,kdens,kx,ky,kz);
kdens_grid(isnan(kdens_grid)) = 0; kdens_grid(isinf(kdens_grid)) = 1;
if SHOWPIX>=2 && doDensComp
    fdens = figure;
    subplot(2,1,1)
        imshow(abs([squeeze(kdens_grid(round(end/2),:,:)), ...
                squeeze(kdens_grid(:,round(end/2),:)), ...
                squeeze(kdens_grid(:,:,round(end/2)))]));
        title("3D grid")
    subplot(2,1,2)
        plot(kdens(1:ndat-2*Nramp));
        title('First platter')
    sgtitle("K-space density compensation")
    drawnow
end

if doT2filter
    if SHOWPIX>=3, fT2 = figure; else, fT2 = []; end
    T2fit = T2filter_220412(T2type,raw,scaninfo,SHOWPIX>=3,fT2);
else
    T2fit = zeros([nframes ndat nleaves nslices ncoils]);
end

% Set up figures and parallel pool options:
fortho = figure; fecho = figure; fcoil = figure; fPD = figure;
if SHOWPIX>=3 || feature('numcores')<=2
    % Run parfor as for (with 0 workers)
    numworkers = 0;
else
    % Use all but 2 to run parfor
    close(fecho,fcoil,fPD)
    numworkers = feature('numcores')-2;
end
if SHOWPIX<1, close(fortho), end
if SHOWPIX >=3 && ~doPhaseDetrend, close(fPD), end

%% Section F: Reconning in loop:
for f = 1:MaxFrame
    
    im_all = zeros(dim3);
    gkdata_all = zeros(dim3);
    
    coilbuf_im = zeros([dim3 ncoils]);
    coilbuf_gk = zeros([dim3 ncoils]);
   
    parfor (c = 1:ncoils, numworkers)
        fprintf('\nReconning frame %d , coil %d\n', f, c);
         
        b = reshape(raw(f,:,c),ndat,nleaves,nslices);
        
        % Do phase detrending if specified
        if doPhaseDetrend && (ncenter > 1)
            b = phaseDetrend_220412(b,ncenter,SHOWPIX>=3,fPD);
            if SHOWPIX>=3
                title(sprintf("Phase Detrend for Frame %d/%d, Coil %d/%d",...
                    f,nframes,c,ncoils))
                drawnow
            end
        end
        
        % Do phase cycling if specified:
        if doPhaseCycle
            if doRots % SERIOS readout
                for ss = 1:2:nslices
                    b(:,:,ss) =  exp(-sqrt(-1)*pi)*(b(:,:,ss))  ;
                end
            else  % stack of spirals:  different order
                for ss = 1:nslices/2
                    b(:,:,ss) =  exp(-sqrt(-1)*pi)*(b(:,:,ss))  ;
                end
            end
        end
        
        % Do T2 filter if specified:
        if doT2filter
            echo_corr = b./squeeze(T2fit(f,:,:,:,c));
            for leaf = 1:nleaves
                b(:,leaf,:) = echo_corr(:,leaf,:)*norm(squeeze(b(:,leaf,:)),1)/norm(squeeze(echo_corr(:,leaf,:)),1);
            end
        end
        
        % Compensating for sampling delay:
        bin  = circshift(b, [Ndel 0 0]); bout = circshift(b, [Ndel2 0 0]);
        bin = bin(1:end/2, :,:);  bout = bout(end/2+1:end, :,:);
        for cnt = 1:size(b,2) % fix the ends
            bin(1:abs(Ndel), cnt) = mean(bin(10:15,cnt),1);
            bout(end-abs(Ndel):end, cnt) = mean(bout(end-10:end-15,cnt),1);
        end
        b = cat(1,bin, bout); % recombine the spirals
        
        % Throwing out ramp/center points
        b = b(Nramp+1:end-Nramp,:,:);
        b(round(end/2)-round(ncenter/2):round(end/2)+round(ncenter/2), :,:) = [];
        
        % Rearranging the order of the echoes to match the acquisition order
        tmpk = [];
        for lf=1 : nleaves
            for sl=1 : nslices
                tmpk = [tmpk ; b(:,lf,sl)];
            end
        end
        tmpk(isnan(tmpk)) = 0;
        
        if SHOWPIX>=3
            figure(fecho)
            plot((abs(tmpk(:)) / (max(abs(tmpk(:))))));
            title(sprintf("Frame %d/%d, Coil %d/%d\nReordered echo train",f,nframes,c,ncoils))
            xlim([1 length(tmpk(:))])
            xlabel('Data points')
            drawnow
        end
        
        % interpolate kspace data along grid:
        gkdata = spiralGrid(ks, tmpk(:), kx, ky, kz);
        gkdata(isnan(gkdata)) = 0; gkdata(isinf(gkdata)) = 0;
        coilbuf_gk(:,:,:,c) = gkdata.*kdens_grid;
        
        % Recon the data:
        if ~doNUFFT % recon using ifft of interpolated grid
            im = DCs(f,c) + fft3d(gkdata.*kdens_grid);
        else % recon using 3D NUFFT
            im = DCs(f,c) + nufft3d(-ks, tmpk(:), fov3, dim3, kdens);
        end
        
        % Store result in a buffer for parallel computing
        coilbuf_im(:,:,:,c) = im;
        
        if SHOWPIX>=3
            figure(fcoil);
            subplot(2,3,1:3)
                imshow(abs([squeeze(im(round(end/2),:,:))', ...
                    squeeze(im(:,round(end/2),:))', ...
                    squeeze(im(:,:,round(end/2)))']));
                axis xy
                caxis('auto'), colormap('gray')
                title("Reconned Image")
            subplot(2,3,4:6)
                imshow(abs(log([squeeze(gkdata(round(end/2),:,:))', ...
                    squeeze(gkdata(:,round(end/2),:))', ...
                    squeeze(gkdata(:,:,round(end/2)))'])));
                axis xy
                caxis('auto')
                title("Gridded k-space data")
            sgtitle(sprintf("Frame %d/%d, Coil %d/%d",f,nframes,c,ncoils))
            drawnow
        end
        
    end
    
    % Combining coil data:
    if ncoils==1
        im_all = coilbuf_im(:,:,:,1);
        gkdata_all = coilbuf_gk(:,:,:,1);
    else
        for c=1:ncoils
            im_all = abs(im_all) + sqrt(-1).*abs(coilbuf_im(:,:,:,c));
            gkdata_all = abs(gkdata_all) + sqrt(-1).*abs(coilbuf_gk(:,:,:,c));
        end
    end
    
    % Flipping image along z-axis:
    im_all = im_all(end:-1:1);
    im_all = reshape(im_all, dim3);
    
    fprintf('\nFrame %d completed.\n' , f);

    if SHOWPIX>=1
        figure(fortho);
        subplot(2,3,1:3)
            imshow(abs([squeeze(im_all(round(end/2),:,:))', ...
                squeeze(im_all(:,round(end/2),:))', ...
                squeeze(im_all(:,:,round(end/2)))']));
            axis xy
            colorbar
            caxis('auto'), colormap('gray')
            title("Reconned Image")
        subplot(2,3,4:6)
            imshow(abs(log([squeeze(gkdata_all(round(end/2)',:,:)), ...
                squeeze(gkdata_all(:,round(end/2),:))', ...
                squeeze(gkdata_all(:,:,round(end/2)))'])));
            axis xy
            caxis('auto')
            title("Gridded k-space data")
        sgtitle(sprintf("Frame %d/%d, All coils combined",f,nframes))
        drawnow
    end
    
    imseries(:,f) =  im_all(:);
    gkseries(:,f) = gkdata_all(:);
    
    imseries4D(:,:,:,f) = im_all;
    gkseries4D(:,:,:,f) = gkdata_all;
    
end

%% Section G: Recon PSF:

% Grid "echo train" containing all ones:
gkdata_psf = spiralGrid(ks, ones(size(ks,1),1), kx, ky, kz);

% Recon the data:
if ~doNUFFT % recon using ifft of interpolated grid
    im_psf = fft3d(gkdata_psf.*kdens_grid);
else % recon using 3D NUFFT
    im_psf = nufft3d(-ks, ones(size(ks,1),1), fov3, dim3, kdens);
end

if SHOWPIX>=2
    
    fpsf = figure;
    subplot(2,3,1:3)
        imshow(abs([squeeze(im_psf(round(end/2),:,:))', ...
            squeeze(im_psf(:,round(end/2),:))', ...
            squeeze(im_psf(:,:,round(end/2)))']));
        axis xy
        caxis('auto'), colormap('gray')
        title("Reconned Image")
    subplot(2,3,4:6)
        imshow(abs(log([squeeze(gkdata_psf(round(end/2),:,:))', ...
            squeeze(gkdata_psf(:,round(end/2),:))', ...
            squeeze(gkdata_psf(:,:,round(end/2)))'])));
        axis xy
        caxis('auto')
        title("Gridded k-space data")
    sgtitle("PSF")
    drawnow
        
end

%% Section H: Write images to nifti file:
if strcmpi(scaleFactor,'auto')
    scaleFactor = DAQfac/max(abs(imseries),[],'all');
end

imseries = scaleFactor * imseries;

write_nii('./kspace_mag.nii', 1e3*log(abs(gkseries(:))), h,0);
write_nii('./timeseries_mag.nii', abs(imseries(:)), h,0);
write_nii('./timeseries_angle.nii', angle(imseries(:))*1000, h,0);
h_psf = h; h_psf.dim(5) = 1;
write_nii('./psf_mag.nii', abs(im_psf(:)),h_psf,0);
write_nii('./psf_ang.nii', angle(im_psf(:)),h_psf,0);
if nframes>3
    hsnr = h;
    hsnr.dim(5) = 1;
    write_nii('./timeseries_snr.nii', ...
        mean(imseries(:,3:end), 2)./ ...
        std(imseries(:,3:end),[],2), hsnr,0);
end

return
%% Section K: Additional functions

function args = spiral3drecon_inputparser(input,defaults)
% INPUTPARSER: function to sort through varargin inputs
    
    % Initilialize output structure to have same format as defaults:
    args = defaults;
    
    % Input validator functions:
    islogic = @(x) (islogical(x) || x==1 || x==0);
    validate = args;
        validate.myPfile = @(x)ischar(x) || @(x)strcmpi(x,'auto');
        validate.saveArgs = @(x)islogic(x);
        validate.zoomFactor = @(x)isnumeric(x);
        validate.scaleFactor = @(x)isnumeric(x) || @(x)strcmp(x,'auto');
        validate.addDC = @(x)islogic(x);
        validate.doNUFFT = @(x)islogic(x);
        validate.doDensComp = @(x)islogic(x) || @(x)strcmpi(x,'auto');
        validate.dcfN = @(x)isnumeric(x);
        validate.dcfR = @(x)isnumeric(x);
        validate.Ndel = @(x)isnumeric(x);
        validate.doRots = @(x)islogic(x) || @(x)strcmpi(x,'auto');
        validate.useKrot = @(x)islogic(x) || @(x)strcmpi(x,'auto');
        validate.doYrot = @(x)islogic(x);
        validate.rotAngle = @(x)isnumeric(x);
        validate.SHOWPIX = @(x)isnumeric(x);
        validate.MaxFrame = @(x)isnumeric(x) || @(x)strcmpi(x,'auto');
        validate.doPhaseDetrend = @(x)islogic(x);
        validate.PDtype = @(x)ischar(x);
        validate.PDspec = @(x)isnumeric(x);
        validate.PDunwrap = @(x)islogic(x);
        validate.doPhaseCycle = @(x)islogic(x);
        validate.doT2filter = @(x)islogic(x);
        validate.T2type = @(x)isnumeric(x);
    
    % Create input parser object:
    p = inputParser;
    
    parmnames = fieldnames(args);
    for i = 1:size(parmnames,1)
        parmname = char(parmnames{i});
        if ~isa(validate.(parmname),'function_handle')
            % In case user forgot to define new parameter validator:
            validate.(parmname) = @(x)1;
        end
        p.addParameter(parmname,defaults.(parmname),validate.(parmname));
    end
    
    p.parse(input{:});
    args = p.Results;
return

function infoarray = spiral3drecon_structinfo(args)
% STRUCTINFO: function to convert information from structure into string
%   array

    fn_args = convertCharsToStrings(fieldnames(args));
    infoarray = [];

    for i = 1:length(fn_args)
        if ischar(args.(fn_args(i)))
            infoarray = [infoarray; sprintf("%s: %s", fn_args(i), args.(fn_args(i)))];
        elseif isnumeric(args.(fn_args(i)))
            infoarray = [infoarray; sprintf("%s: %d", fn_args(i), 1*args.(fn_args(i)))];
        elseif isstruct(args.(fn_args(i)))
            fn_fieldstruct = convertCharsToStrings(fieldnames(args.(fn_args(i))));
            for j = 1:length(fn_fieldstruct)
                if ischar(args.(fn_args(i)).(fn_fieldstruct(j)))
                    infoarray = [infoarray; sprintf("%s.%s: %s", fn_args(i), fn_fieldstruct(j), ...
                        args.(fn_args(i)).(fn_fieldstruct(j)))];
                elseif isnumeric(args.(fn_args(i)).(fn_fieldstruct(j)))
                    infoarray = [infoarray; sprintf("%s.%s: %.2f", fn_args(i), fn_fieldstruct(j), ...
                        1*args.(fn_args(i)).(fn_fieldstruct(j)))];
                end
            end
        end
    end
    
return
