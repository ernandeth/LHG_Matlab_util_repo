function [im_4D,ksp_4D,psf_3D,args] = recon_vsasl3dflex(varargin)
% recon_vsasl3dflex reconstruct image data generated from vsasl3dflex
% sequence by loading k-space trajectory and data, then reconstructing
% on a frame/coil-wise basis using the non-uniform FFT algorithm (Jeff
% Fessler - MIRT)
%
%   Luis Hernandez-Garcia, David Frey @UM 2022
%
% recon_vsasl3dflex() will reconstruct data from working directory using
% default parameters. Several files are required for this to run:
%
%   - Pfile (i.e. 'P*.7')
%   - First platter gradient data (i.e. 'grad.txt')
%   - View parameter tables (i.e. 'kviews.txt')
%   - Spherical trajectory data (i.e. 'ktraj_sph.txt')
%
% Several parameters can be passed using pair-wise varargin formatting
% (i.e. recon_vsasl3dflex('<arg name>',<arg value>) ). If recon has been
% ran previously in same working directory, the arguments may be saved t o
% the file 'lastreconargs.mat', and can be reused by calling 'lastargs' as
% first argument
% (i.e. recon_vsasl3dflex('lastargs','<arg name>',<arg value>) ). Other
% input parameters include:
%   
%   Basic File Info:
%       myPfile
%           - name of pfile (str/char) to load for reconstruction
%           - if 'auto' is passed, recon will use first Pfile in working
%               directory
%           - default is 'auto'
%       raw
%           - raw dataset (<nframes> x <total datapts> x <ncoils> complex)
%           - if 'auto' is passed, recon will grab raw data from Pfile.
%           - data can be manually passed if there is some preprocessing
%               the user wants to do (despiking, etc.)
%           - if raw is manually passed, user must also manually pass the
%               scaninfo structure
%           - note: if you want to despike the data prior to passing into
%               recon, the recon will automatically checked for data saved
%               to a file
%           - default is 'auto'
%       scaninfo
%           - scanner information structure (struct)
%           - if 'auto' is passed, recon will grab scaninfo from Pfile.
%           - must be manually passed if raw is also manually passed
%           - default is 'auto'
%
%   Debug setting:
%       SHOWPIX
%           - debugging level (int) to show figures as image is
%               reconstructed
%           - a couple different levels can be passed
%               (0): Show no figures, run as normal and write to files at 
%                   the end
%               (>=1): Show orthogonal views of frame by frame
%                   reconstruction
%               (>=2): Show raw echo train for example leaf/frame/coil,
%                   k-space trajectory (vector components and 3D plot), and
%                   orthogonal views & curve of density compensation
%                   function
%               (>=3): Show simulated view of readout as expected to be
%                   displayed on an oscilloscope, orthogonal views of point
%                   spread function, and temporal SNR map if nframes>=3
%               (>=4): Show orthogonal views of coil by coil
%                   reconstruction,
%                   <<<<<<<<<< to finish later >>>>>>>>>>>>>>>>
%                   WARNING: running with this setting requires that
%                   reconstruction loop be ran with 0 workers to display
%                   figures, meaning it will be much slower
%           - default is 0
%       maxFrame
%           - maximum number of frames to recon (int)
%           - if 'auto' is passed, all frames will be reconned
%           - default is 'auto'
%       numWorkers
%           - number of cores (workers) to use for parallel pool (int)
%           - default is total number of cores available
%       saveArgs
%           - option to save args structure to file (bool)
%           - default is 1
%       saveKspace
%           - option to grid and save kspace data (bool)
%           - takes significantly more time
%           - default is 0
%
%   Basic Image Properties:
%       resFactor
%           - resolution interpolated upsampling factor (double)
%           - (i.e. resFactor = 2 will result in 2x the output resolution)
%           - default is 1
%       zoomFactor
%           - zooming factor (double)
%           - (i.e. zoomFactor = 2 will result in an image zoomed in by 2x)
%           - default is 1
%       doScl
%           - option to save full dynamic range images to nifti (bool)
%           - scales images from -32767 to 32767 before saving to .nii and
%               saves scaling parameters scl_slope & scl_inter to header
%           - default is 1
%       M0frames
%           - number of M0 frames (int)
%           - used for frame outlier removal and tsnr calculation
%           - default is 2
%
%   Reconstruction Tuning Parameters:
%       doNUFFT
%           - option to use NUFFT to reconstruct kspace data (bool)
%           - default is 1
%       Ndel
%           - sampling delay correction (int)
%           - default is -1
%       Nramp
%           - number of ramp points to throw out (int)
%           - if 'auto' is passed, recon will determine based on
%               radial kspace data from ktraj_sph.txt
%           - default is 'auto'
%       Nechoes
%           - override for number of echoes (slices) to reconstruct in each
%               train (int)
%           - if 'auto' is passed, recon will use all echoes
%           - default is 'auto'
%       orderPhaseDetrend
%           - option to correct off-resonance effects by fitting a
%               polynomial to phase drift at navigator points, the value
%               specifies the order of the fit (int)
%           - if -1 is passed, no phase detrending will be performed
%           - default is -1
%       T2weight
%           - weighting factor for T2 correction (float)
%           - a lower weight will result in a less severe T2 correction by
%               multiplying the fitted rate
%           - default is 0 (no correction)
%
%       Density Compensation Function (dcf) Options:
%       dcfRecycle
%           - option to save/load dcf from a previous iteration if
%               parameters are consistent (bool)
%           -default is 1
%       dcfP
%           - power of exponential distance weighting function (double)
%           - i.e. density is computed as the mean(<distances to dcfN
%               closest neighbors>)^dcfP
%           - default is 2
%       dcfN
%           - number of neighbors to use for density calculation (int)
%           - refer to calculation of density in documentation for dcfP
%               directly above
%           - default is 50
%       despike
%           - option to despike k-space data along frame dimension (bool)
%           - recommend to despike the data and pass in manually before
%               running the recon (using despike_vsasl3dflex.m)
%           - default is 0
%       despikeRecycle
%           - option to save/load despiked data from a previous iteration
%               if parameters are consistent
%       rf_phase_cycle
%           - option to correct for rf phase cycling (90,180,-180,180,-180...)
%           - default is 0
%
%   Scanner Constants:
%       gamma
%           - gyromagnetic ratio (float)
%           - default is ~4257 Hz/s/Gauss
%       DACFactor
%           - A/D conversion ratio (int)
%           - default is 32767 (or 2^15 - 1)
%
% Several output parameters can be returned by
% [args,ksp_4D,im_4D,psf_3D] = recon_vsasl3dflex(); These include (in
% order):
%
%       im_4D
%           - 4D array containing all reconstructed image data (4D complex
%               double)
%       ksp_4D
%           - 4D array containing all k-space data (4D complex double)
%       psf_3D
%           - 3D array containing reconstructed point spread function (3D
%               complex double)
%       args
%           - structure containing all arguments passed (struct)
%
% Several output files can be generated upon reconstruction, including:
%
%   Nifti Image Files
%       timeseries_mag.nii
%           - magnitude image of reconstructed timeseries
%       timeseries_ang.nii
%           - phase image of reconstructed timeseries
%       kspace_mag.nii
%           - magnitude image of k-space timeseries data
%       kspace_ang.nii
%           - phase image of k-space timeseries data
%       psf_mag.nii
%           - magnitude image of point spread function
%       psf_ang.nii
%           - phase image of point spread function
%       tSNR.nii
%           - magnitude image of temporal SNR map
%           - will only be generated when more than 3 frames are
%               reconstructed
%
%   Documentation Files
%       history.txt
%           - log of all parameters and dates of past recons for current
%               working directory
%       lastreconargs.mat
%           - .mat file containing used arguments for future
%               reproducibility (see explanation of 'lastargs' above)
%           - if 'lastargs' is used as first input argument, contents of
%               this file will be used as defaults
%
%

% Start timer
tStart = tic;

%% SECTION A: Defining variables and setting defaults
% Note: To add a parameter, add them into the structure below along with
%   its default value. Then add the line
%   "<name of parm> = args.<name of parm>" under "Convert fields from args
%   into variable" in section B.

    % Set defaults for variable input parameters in structure below:
    defaults = struct(...
    ... % Basic file info:
        'myPfile', 'auto', ... % 'auto' sets name of pfile based on first in working directory
        'raw', 'auto', ... % 'auto' grabs raw data from pfile
        'scaninfo', 'auto', ... % 'auto' grabs scaninfo from pfile
    ...
    ... % Debug settings:
        'SHOWPIX', 0, ...
        'maxFrame', 'auto', ...
        'numWorkers', feature('numcores'), ...
        'saveArgs', 1, ...
        'saveKspace', 0, ...
    ...
    ... % Basic image properties:
        'resFactor', 1, ...
        'zoomFactor', 1, ...
        'doScl', 1, ...
        'M0frames', 2, ...
    ...
    ... % Reconstruction tuning parameters:
        'doNUFFT', 1, ...
        'Ndel', -1, ...
        'Nramp', 'auto', ...
        'Nechoes', 'auto', ...
        'orderPhaseDetrend', -1, ...
        'T2weight', 0, ...
        'echoLowpass', 1, ...
        'dcfRecycle', 1, ...
        'dcfP', 2, ...
        'dcfN', 50, ...
        'despike', 0, ...
        'despikeRecycle', 1, ...
        'rf_phase_cycle', 0, ...
    ...
    ... % Scanner constants:
        'gamma', 267.5222 * 1e6 * 1e-4 / 2 / pi, ...
        'DACFactor', 2^15-1 ...
        );
    
%% SECTION B: Parsing variable input parameters

    fprintf('\nParsing input parameters...');

    % If user calls 'lastargs' as first variable input, replace defaults
    % with last saved args
    if numel(varargin)>0 && strcmpi(varargin{1},'lastargs')
        lastreconargs = load('lastreconargs.mat','args');
        defaults = lastreconargs.args;
        varargin(1) = [];
    end
    
    % Parse through variable inputs using matlab's built-in input parser
    p = inputParser;    
    parmnames = fieldnames(defaults);
    for i = 1:size(parmnames,1)
        parmname = char(parmnames{i});
        p.addParameter(parmname,defaults.(parmname),@(x)1);
    end
    p.parse(varargin{:});
    args = p.Results;

    % If variables are defined with auto flag, handle them using function
    isauto = @(x)strcmpi(x,'auto');
    
    % Auto assign Pfile
    if isauto(args.myPfile)
        % Assign pfile based on name of first pfile in directory
        pfiles = dir('P*.7');
        args.myPfile = pfiles(1).name;
    end
    
    % Convert fields from args into variables
    myPfile 	= args.myPfile;
    raw 	= args.raw;
    scaninfo 	= args.scaninfo;
    SHOWPIX 	= args.SHOWPIX;
    maxFrame 	= args.maxFrame;
    numWorkers 	= args.numWorkers;
    saveArgs 	= args.saveArgs;
    saveKspace 	= args.saveKspace;
    resFactor 	= args.resFactor;
    zoomFactor 	= args.zoomFactor;
    doScl 	= args.doScl;
    M0frames 	= args.M0frames;
    doNUFFT 	= args.doNUFFT;
    Ndel 	= args.Ndel;
    Nramp 	= args.Nramp;
    Nechoes 	= args.Nechoes;
    orderPhaseDetrend = args.orderPhaseDetrend;
    T2weight 	= args.T2weight;
    echoLowpass = args.echoLowpass;
    dcfRecycle 	= args.dcfRecycle;
    dcfP 	= args.dcfP;
    dcfN 	= args.dcfN;
    despike 	= args.despike;
    despikeRecycle = args.despikeRecycle;
    rf_phase_cycle = args.rf_phase_cycle;
    gamma 	= args.gamma;
    DACFactor 	= args.DACFactor;

    % Validate number of workers for parfor loop
    if numWorkers > feature('numCores')
        warning('\nSpecified number of workers (%d) to use for parfor exceeds number of cores (%d)', ...
            numWorkers,feature('numCores'));
        numWorkers = feature('numCores');
    elseif numWorkers < 0
        warning('\nSpecified number of workers (%d) cannot be negative', ...
            numWorkers);
        numWorkers = 0;
    end
    
    % Start parallel pool
    if numWorkers>0 && isempty(gcp('nocreate'))
        fprintf('\n');
        parpool(numWorkers);
    end
    
%% SECTION C: Extracting and configuring aquisition information

    fprintf('\nExtracting/Configuring acquisition information...');

    % Load data from files
    if isauto(raw)
        [~,raw,scaninfo] = evalc('read_raw_3d(myPfile,0);'); % Raw data (Pfile)
    elseif isauto(scaninfo)
        error('if raw is manually passed, you must also pass scaninfo structure');
    elseif strcmpi(raw,'manual')
        warning('Last recon manually passed raw, must manually pass again if you want to use the same raw');
        fprintf('\nProceeding with reading raw from a pfile...');
        [~,raw,scaninfo] = evalc('read_raw_3d(myPfile,0);'); % Raw data (Pfile)
    else
        args.raw = 'manual';
        args.scaninfo = 'manual';
    end
    ktraj_sph = load('ktraj_sph.txt'); % Spherical coordinates of k-space traj.
    g = load('grad.txt'); % Gradient waveform for first platter
    kviews = load('kviews.txt'); % View parameter tables
    
    % Determine if RO is a stack of spirals based on rotation angles
    if (all(kviews(:,4)==0) && all(kviews(:,5)==0))
        isSOS = 1;
    else
        isSOS = 0;
    end
    
    % Get info from scaninfo
    dt = scaninfo.ts;
    nslices = scaninfo.nslices; % number of slices
    nleaves = scaninfo.npr; % number of interleaves
    nframes = scaninfo.nphases; % number of frames
    ndat = scaninfo.ndat; % number of data points per echo
    ncoils = scaninfo.ncoils; % number of coils used
    slthick = scaninfo.slthick/10; % slice thickness (cm)
    te = scaninfo.te; % echo time (ms)
    tr = scaninfo.tr*1e-3; % scan repetition time (s)
    fov = ones(1,3)*scaninfo.opfov/zoomFactor; % field of view (cm)
    dim = ones(1,3)*ceil(resFactor*scaninfo.opxres); % image dimension (number of voxels)
    if isSOS
        fov(3) = slthick*nslices/zoomFactor;%%%%%%%%%%%%%%%%%%%
        dim(3) = ceil(resFactor*nslices*nleaves);
    end
    
    % Auto assign maxFrame
    if isauto(maxFrame) && SHOWPIX>=4
        args.maxFrame = 1; maxFrame = args.maxFrame;
    elseif isauto(maxFrame)
        args.maxFrame = nframes; maxFrame = args.maxFrame;
    end
    
    % Perform frame-wise k-space despiking
    if despike
        raw_corr = despike_vsasl3dflex(raw,M0frames,despikeRecycle);
        raw = raw_corr;
    end
    
    % Determine Kmax and create a cartesian grid:
    Kmax = dim./fov/2;
    [Kx,Ky,Kz]= meshgrid( ...
        linspace(-Kmax(1),Kmax(1),dim(1)), ...
        linspace(-Kmax(2),Kmax(2),dim(2)), ...
        linspace(-Kmax(3),Kmax(3),dim(3)));
    
    % Reshape raw into individual echoes
    raw = reshape(raw,[nframes,ndat,nleaves,nslices,ncoils]);
    
    % Attenuate unwanted echoes and change nslices
    if isauto(Nechoes) || ~(Nechoes>0 && Nechoes<=nslices)
        args.Nechoes = nslices;
        Nechoes = args.Nechoes;
    end
    nslices = Nechoes;
    raw = raw(:,:,:,1:nslices,:);
    raw_filt = raw;
    
    % Perform Lowpass filtering
    if echoLowpass < 1
        raw_filt = lowpass_vsasl3dflex(raw,echoLowpass);
    end
    
    % Plot raw echo train at middle of sequence
    if args.SHOWPIX>=2 
        genfig_rawecho(raw,raw_filt,echoLowpass);
    end
    
    raw = raw_filt;
    
%% SECTION D: Processing k-space trajectory/density compensation

    fprintf('\nProcessing k-space trajectory...');

    if isauto(Nramp)
        % Find maxima in r curve of k-space trajectory
        [~,Nramp] = min(abs(ktraj_sph(1:round(end/4),1) - Kmax(1)));
    end
        
    % Determine navigator points from spherical trajectory
    % (center of k-space should be consistent from echo to echo)
    navpts = find(ktraj_sph(:,1) == 0);
    % ... and reject ramp points that may be at zero kspace
    navpts(navpts>ndat-Nramp) = []; 
    navpts(navpts<Nramp) = [];
    ncenter = length(navpts);
    
    % Check that acquisition size is equal to trajectory size
    if (ndat<size(g,1))
        fprintf('\n');
        warning(['acquistion length (ndat = %d) is less', ...
            ' than gradient length (size(g,1) = %d)\n\tClipping size(g,1) to ndat'], ...
            ndat, size(g,1));
        g = g(1:ndat,:);
    elseif (ndat>size(g,1))
        fprintf('\n');
        warning(['acquistion length (ndat = %d) is greater', ...
            ' than gradient length (size(g,1) = %d)\n\tClipping ndat to size(g,1)'], ...
            ndat, size(g,1));
        ndat = size(g,1);
        raw = raw(:,1:ndat,:,:,:);
    end
    
    % Calculate kspace trajectory for first spiral
    ks = zeros(ndat,3,nleaves,nslices);
    ks_0 = gamma*dt*(cumsum(g));
    
    % Translate first spiral to list of views from kviews.txt  
    for leafn = 1:nleaves
        for slicen = 1:nslices
            
            % Calculate overall view (spiral) index
            spiraln = (leafn-1)*nslices + slicen;
            
            % Transform trajectory by applying transformation matrix
            T = reshape(kviews(spiraln,end-8:end),3,3)';
            
            ks(:,:,leafn,slicen) = (T*ks_0')';
            
        end
    end
    
    % Plot k-space trajectory
    if SHOWPIX>=2
        genfig_ktraj(ks);
    end
    
    % Plot gradient/RF scope view with corresponding echo train
    if SHOWPIX>=3
        % ISSUE TO FIX:
        %genfig_ROscope(ks,raw,TE,dt);
    end
    
%% SECTION E: Pre-recon loop

    fprintf('\nSetting up for recon loop...');
    
    % Initialize output arrays
    ksp_4D = zeros([dim maxFrame]);
    im_4D = zeros([dim maxFrame]);
    
    % Store DC values (data at center of kspace)
    DCs = squeeze(mean(raw(:,navpts,1,1,:),2));
    
    % Perform phase detrending
    if (orderPhaseDetrend > -1)
        fprintf('\n\tPerforming phase detrending...');
        [raw_corr,PDfits] = phaseDetrend_vsasl3dflex(raw,navpts,orderPhaseDetrend,isSOS);
        % Show phase detrending
        if SHOWPIX>=3
            genfig_PD(raw,raw_corr,PDfits,orderPhaseDetrend)
        end
        raw = raw_corr;
    end
    
    % Perform T2 filtering
    if T2weight>0
        fprintf('\n\tPerforming T2 filtering with weight %.2f...',T2weight)
        [raw_corr,T2,T2curve] = T2filter_vsasl3dflex(raw,T2weight,navpts,te,dt);
        % Show T2 filtering
        if SHOWPIX>=3
            genfig_T2(raw,raw_corr,T2,T2curve)
        end
        raw = raw_corr;
    end
    
    % Correct for 180 phase cycling
    if rf_phase_cycle
        raw = raw .* reshape(exp((-1).^(1:nslices)*sqrt(-1)*pi/2),[1 1 1 nslices 1]);
    end
    
    % Correct sampling delay
    if abs(Ndel)>0
        fprintf('\n\tCorrecting sampling delay by %d samples...',Ndel);
        raw_tmp = circshift(raw, [0 Ndel 0 0 0]);
        raw_in = raw_tmp(:,1:round(end/2),:,:,:);
        raw_out = raw_tmp(:,round(end/2)+1:end,:,:,:);
        % ... and fix the ends of the spiral
        raw_in(:,1:abs(Ndel),:,:,:) = repmat(mean(raw_in(:,abs(Ndel):abs(Ndel)+5,:,:,:),2),[1 abs(Ndel) 1 1 1]);
        raw_out(:,end-abs(Ndel)+1:end,:,:,:) = repmat(mean(raw_out(:,end-abs(Ndel)-5:end-abs(Ndel),:,:,:),2),[1 abs(Ndel) 1 1 1]);
        raw = cat(2,raw_in,raw_out);
    end
    
    % Remove ramp points
    if Nramp>0
        fprintf('\n\tRemoving %d ramp points...', Nramp);
        raw(:,1:Nramp,:,:,:) = []; raw(:,end-Nramp:end,:,:,:) = [];
        ks(1:Nramp,:,:,:) = []; ks(end-Nramp:end,:,:,:) = [];
    end
    
    % Compute density compensation function
    if doNUFFT
        dcf = DCF_vsasl3dflex(ks,dcfP,dcfN,dcfRecycle,numWorkers);
    elseif numWorkers > 0
        dcf = []; % Need to at least declare it as a variable for parfor loop
    end
    
    % Show density compensation function
    %{
    if SHOWPIX>=2 && ~doNUFFT
        genfig_dcf(ks,Kx,Ky,Kz,dcf);
    end
    %}
    % Extract individual trajectory components for gridding
    ksx = reshape(ks(:,1,:,:),[],1);
    ksy = reshape(ks(:,2,:,:),[],1);
    ksz = reshape(ks(:,3,:,:),[],1);
    
    % if SOS, dims of neighborhood should be 6x6x1
  %{
         nufft_args = {dim,... % size of output
        %10*ones(3,1),...
        [6 6 1]',...        % size of neighborhood for interpolation. eg-  6 6 6 works well.  
                            % BUT note that if already cartesian, you don't need a neighborhood
        %2*dim,...
        [2 2 1].*dim, ...   % oversampling factor for FFT.  you don't need to oversample if it's already cartesian
        dim/2,...           % this determines where 0 pixel is - thinkg fftshift
        'table',...         % how to interpolate:  accuracy v. memory. 
                            %   'table' precomputes a table of kernel weights
        2^10,...            %  how many samples in the table.
        'minmax:kb'};       % kernel choice 
  %}
    %%
    % in SOS we expect the dcf to be the same in all platters
    Nnbrs = 6*ones(3,1);
    oversamp = 2*dim;
    
    if isSOS
       Nnbrs(3) = 1;
       oversamp=[2 2 1].*dim;
    end
    % Set up nufft
    nufft_args = {dim,...
        Nnbrs,...
        oversamp, ...
        dim/2,...   
        'table',... 
        2^10,...   
        'minmax:kb'}; 
    omega = -[ksx ksy ksz] * 2*pi .* fov ./ dim;
    Gob = Gnufft(true(dim), [{omega}, nufft_args(:)']);
    
    % Reconstruct the point spread function
    fprintf('\n\tReconstructing point spread function...');
    if doNUFFT
        psf_3D = reshape(Gob' * (dcf(:) .* ones(size(ksx))), dim);
    else
        %{
        warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
        
        F = scatteredInterpolant(-ksx,-ksy,-ksz,ones(size(ksx)),'linear','none');
        ksp = F(Kx,Ky,Kz);
        ksp(isnan(ksp)) = 0; ksp(isinf(ksp)) = 0;
        psf_3D = fft3d(ksp);
        %}
        fprintf('(using grid3d_lhg)')
        [psf_3D dcf kernel] = grid3d_lhg(...
            ksx, ksy, ksz, ...
            ones(size(ksx)), ...
            Kx, Ky, Kz, ...
            2, dcf);
    end
    
    % Plot psf
    if SHOWPIX>=3
        genfig_psf(psf_3D,fov)
    end
    
%% SECTION F: Recon loop

    fprintf('\nBeginning recon loop...');
    
    % Warn user of intentionally slow recon loop due to debug level
    if SHOWPIX>3
        warning('\nRecon loop running with zero workers since SHOWPIX = %d > 3',...
            SHOWPIX);
    end
    
    for framen = 1:maxFrame
        fprintf('\nReconning frame %d/%d, ', framen,maxFrame);
        if usejava('desktop')
            progbar = sprintf('\b|\n');
            fprintf(['Coil-wise progress:\n' repmat('.',1,ncoils) ' (%d coils total)\n\t'], ncoils);
        else
            progbar = '';
            fprintf('(cannot display coil-wise progress in nodesktop mode) ');
        end
        
        % Set coil-wise data buffers to zeros
        im_coilbuf = zeros([dim ncoils]);
        ksp_coilbuf = zeros([dim ncoils]);
        
        parfor (coiln = 1:ncoils, numWorkers*(SHOWPIX<4))
            
            % Get data for individual coil/frame
            data = raw(framen,:,:,:,coiln);
            
            % Grid the k-space data using simple linear interpolation
            if (saveKspace || ~doNUFFT || SHOWPIX>=0)
                %{
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
                F = scatteredInterpolant(ksx,ksy,ksz,data(:),'nearest','none');
                ksp = F(Kx,Ky,Kz); %#ok<PFTIN>
                ksp(isnan(ksp)) = 0; ksp(isinf(ksp)) = 0;
                % ... and store it into the coil-wise data buffer
                ksp_coilbuf(:,:,:,coiln) = ksp;
                %}
                        
                fprintf('using grid3d_lhg...coil %d ', coiln)

                [ksp dcf kernel] = grid3d_lhg(...
                    ksx, ksy, ksz, ...
                    data(:), ...
                    Kx, Ky, Kz, ...
                    2, []);
                
                ksp_coilbuf(:,:,:,coiln) = ksp;

            end
            
            % Perform non-uniform fast fourier transform
            if doNUFFT
                im = DCs(framen,coiln) + reshape(Gob' * (dcf(:) .* data(:)), dim);
            else
                im = fftshift(fftn(ksp,dim));
            end
            % ... and store it into the coil-wise data buffer
            im_coilbuf(:,:,:,coiln) = im;
            
            % Show reconstruction for current coil
            if SHOWPIX>=4
                genfig_ovcoil(im,ksp,fov,Kmax,framen,nframes,coiln,ncoils);
            end
            
            fprintf('%s',progbar);
        end
        
        % Combine coils
        im = sqrt(mean(abs(im_coilbuf).^2,4));
        ksp = sqrt(mean(abs(ksp_coilbuf).^2,4));
        
        if ~doNUFFT
            fk = fftshift(abs(fftn(kernel)));
            fk = fk- min(fk(:));
            fk = fk/max(fk(:));
            im = im./(0.01+fk);
        end
        
        im_4D(:,:,:,framen) = im;
        ksp_4D(:,:,:,framen) = ksp;
        
        % Show reconstruction for current frame
        if SHOWPIX>=1
            genfig_ovframe(im,ksp,fov,Kmax,framen,nframes);
        end
        
        fprintf('\b (Done!)');
        
    end
    
    % Calculate tSNR map:
    if maxFrame > M0frames
        tSNR_3D = squeeze(mean(im_4D(:,:,:,(M0frames+1):end),4)) ./ ...
            squeeze(std(im_4D(:,:,:,(M0frames+1):end),[],4));
    end
    
    % Show tSNR map:
    if SHOWPIX>=3 && maxFrame > M0frames
        genfig_tSNR(tSNR_3D,fov);
    end
    
%% SECTION G: Configuring and writing to output

    fprintf('\nConfiguring and writing to output...');
    
    % Save nifti files
    if saveKspace
        matrix2nii('./kspace_mag.nii',abs(ksp_4D),1./fov,tr,doScl);
    end
    matrix2nii('./timeseries_mag.nii',abs(im_4D),fov,tr,doScl);
    matrix2nii('./psf_mag.nii',abs(psf_3D),fov,tr,doScl);
    matrix2nii('./psf_ang.nii',angle(psf_3D),fov,tr,doScl);
    if maxFrame > M0frames
        matrix2nii('./tSNR.nii',abs(tSNR_3D),fov,tr,doScl);
    end
    
    if saveArgs
        save lastreconargs.mat args
    end
    
    tEnd = toc(tStart);
    fprintf('\nRecon complete. Elapsed time = %fs\n',tEnd)
    
end



%% APPENDIX A: Mid-level helper functions

function cfigopen(figname)

    % Check if figure with specified name is open
    if isempty(findobj('type','figure','name',figname))
        % If not, open it
        figure('name',figname);
    else
        % If so, make it the current figure
        figure(findobj('type','figure','name',figname))
    end
    
end

function showov(im,rowtitle,rown,nrows,bxlabel,bylabel,bzlabel)

    % Normalize image
    if max(abs(im(:)))>0
        im = im./max(abs(im(:)));
    end

    % Plot cut through center of first dimension
    subplot(nrows,3,(rown-1)*3+1)
        imshow(squeeze(abs(im(round(end/2),end:-1:1,end:-1:1)))')
        axis on
        xlabel('Y'); ylabel('Z');
        xticks([1 size(im,2)]), xtickangle(0)
        xticklabels(bylabel)
        yticks([1 size(im,3)]), ytickangle(90)
        yticklabels(bzlabel(end:-1:1))
        
    % Plot cut through center of second dimension
    subplot(nrows,3,(rown-1)*3+2)
        imshow(squeeze(abs(im(end:-1:1,round(end/2),end:-1:1)))')
        axis on
        xlabel('X'); ylabel('Z');
        xticks([1 size(im,1)]), xtickangle(0)
        xticklabels(bxlabel)
        yticks([1 size(im,3)]), ytickangle(90)
        yticklabels(bzlabel(end:-1:1))
        title(rowtitle)
        
    % Plot cut through center of third dimension
    subplot(nrows,3,(rown-1)*3+3)
        imshow(squeeze(abs(im(end:-1:1,end:-1:1,round(end/2))))')
        axis on
        xlabel('X'); ylabel('Y');
        xticks([1 size(im,1)]), xtickangle(0)
        xticklabels(bxlabel)
        yticks([1 size(im,2)]), ytickangle(90)
        yticklabels(bylabel(end:-1:1))
        
end



%% APPENDIX B: Low-level figure generating functions

function genfig_rawecho(raw,raw_filt,echoLowpass)

    cfigopen('Raw Echo Train')
    
    nframes = size(raw,1);
    ndat = size(raw,2);
    nleaves = size(raw,3);
    nslices = size(raw,4);
    ncoils = size(raw,5);

    subplot(2,1,1)
        plot(...
            reshape(abs(raw(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1));
        ylabel('Magnitude')
        xlim([0 ndat*nslices]);
        xticks(ndat*(1/2 + 0:nslices));
        xticklabels(1:nslices);
        xlabel('Slice (echo) index')
        
        if echoLowpass < 1
            hold on
            plot(...
                reshape(abs(raw_filt(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1));
            hold off
            legend('Unfiltered',sprintf('Lowpass (w_c = %.3f) Filtered',echoLowpass))
        end

    subplot(2,1,2)
        plot(...
            reshape(angle(raw(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1));
        ylabel('Phase')
        yticks([-pi 0 pi]);
        yticklabels(["-\pi" "0" "\pi"]);
        xlim([0 ndat*nslices]);
        xticks(ndat*(1/2 + 0:nslices));
        xticklabels(1:nslices);
        xlabel('Slice (echo) index')

    sgtitle(sprintf('Raw echo train for leaf %d, frame %d, coil %d',...
        ceil(nleaves/2),ceil(nframes/2),ceil(ncoils/2)));

    drawnow
    
end

function genfig_ktraj(ks)
    
    cfigopen('K-space Trajectory')

    ndat = size(ks,1);
    nleaves = size(ks,3);
    nslices = size(ks,4);
    
    subplot(2,3,[1 4])
        x = 1:size(ks,1);
        spacing = max(abs(ks(:)))*2;
        plot(...
            x,ks(:,1,ceil(nleaves/2),ceil(nslices/2))+spacing,...
            x,ks(:,2,ceil(nleaves/2),ceil(nslices/2)),...
            x,ks(:,3,ceil(nleaves/2),ceil(nslices/2))-spacing)
        [minx,maxx] = bounds(x); xlim([minx maxx]), xticks([])
        yticks(spacing*(-1:1)), yticklabels(["K_z" "K_y" "K_x"])
        title(sprintf("Spiral vector components\nfor leaf %d, slice %d",...
            ceil(nleaves/2),ceil(nslices/2)))
        
    subplot(2,3,[2 3 5 6])
        plot3(...
            reshape(ks(:,1,:,:),[],1),...
            reshape(ks(:,2,:,:),[],1),...
            reshape(ks(:,3,:,:),[],1),...
            'LineWidth',0.5)
        hold on
        plot3(...
            reshape(ks(:,1,ceil(nleaves/2),ceil(nslices/2)),[],1),...
            reshape(ks(:,2,ceil(nleaves/2),ceil(nslices/2)),[],1),...
            reshape(ks(:,3,ceil(nleaves/2),ceil(nslices/2)),[],1),...
            'LineWidth',2)
        hold off
        legend("Full trajectory",sprintf("Spiral for leaf %d, slice %d",...
            ceil(nleaves/2),ceil(nslices/2)),"Location","SouthOutside")
        xlabel('K_x'), ylabel('K_y')
        title("3D view")
        
    sgtitle("K-space Trajectory");

    drawnow
    
end

function genfig_ROscope(ks,raw,TE,dt)
    
    cfigopen('Readout Oscilloscope View')
    
    ndat = size(ks,1);
    nleaves = size(ks,3);
    nslices = size(ks,4);

    ndat_all = round(TE*1e-3/dt);
    if (mod(ndat_all-ndat+1,2)~=0), ndat_all = ndat_all+1; end
    
    g = padarray(diff(ks,1,1),[(ndat_all-ndat+1)/2,0,0,0],0);
    
    spacing = max(abs(g(:)))*2;
    
    RF = zeros(ndat_all,nleaves,nslices);
    RF(1,:,2:end) = max(abs(g(:)));
    RF(1,:,1) = max(abs(g(:)))/2;
    
    echo = squeeze(abs(raw(1,:,:,:,1)))...
        *1.5*max(abs(g(:)))/max(abs(raw(1,:,:,:,1)),[],'all');
    echo = padarray(echo(1:end-1,:,:),[(ndat_all-ndat+1)/2,0,0,0],0);
    
    for leafn = 1:nleaves
        
        gxl = reshape(g(:,1,leafn,:),[],1) + 2*spacing;
        gyl = reshape(g(:,2,leafn,:),[],1) + 1*spacing;
        gzl = reshape(g(:,3,leafn,:),[],1) - 0*spacing;
        RFl = reshape(RF(:,leafn,:),[],1) - 1*spacing;
        echol = reshape(echo(:,leafn,:),[],1) - 2*spacing;

        x = (ndat_all*nslices*(leafn-1):ndat_all*nslices*leafn-1) + (leafn-1)*ndat_all;
        plot(x,gxl,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
        hold on
        plot(x,gyl,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
        plot(x,gzl,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
        plot(x,RFl,'Color',[0.4940 0.1840 0.5560],'LineWidth',2);
        plot(x,echol,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
        xline(ndat_all*nslices*(leafn-1) + (leafn-1)*ndat_all,...
            'Label',sprintf('Leaf %d',leafn),'LineStyle','--');
        xline(ndat_all*nslices*leafn-1 + (leafn-1)*ndat_all,...
            'LineStyle','--');
            
    end
    
    hold off
    xlim(ndat_all*[-0.5 (nslices+1)*nleaves+0.5]), xticks([])
    yticks(spacing*(-2:2)), yticklabels(["Echo" "RF" "G_z" "G_y" "G_x"])
    sgtitle("Readout oscilloscope view for frame 1, coil 1");
    
    drawnow
    
end

function genfig_dcf(ks,Kx,Ky,Kz,dcf)
    
    cfigopen('Density Compensation Function')
    
    ndat = size(ks,1);
    nleaves = size(ks,3);
    nslices = size(ks,4);
    
    ksx = reshape(ks(:,1,:,:),[],1);
    ksy = reshape(ks(:,2,:,:),[],1);
    ksz = reshape(ks(:,3,:,:),[],1);
    
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
    dcf_grid = griddata(ksx,ksy,ksz,dcf(:),Kx,Ky,Kz);
    
    showov(dcf_grid, 'K-space', 1, 2, ...
        [sprintf("%.1f",min(ksx)) sprintf("%.1f",max(ksx))], ...
        [sprintf("%.1f",min(ksy)) sprintf("%.1f",max(ksy))], ...
        [sprintf("%.1f",max(ksz)) sprintf("%.1f",min(ksz))]);
        
    subplot(2,3,4:6)
        plot(dcf(:,ceil(nleaves/2),ceil(nslices/2)))
        title(sprintf('Weighting curve for leaf %d, slice %d',ceil(nleaves/2),ceil(nslices/2)))
        xlim([1 ndat])
        
    sgtitle('Density compensation function');

    drawnow
    
end

function genfig_PD(raw,raw_corr,PDfits,orderPhaseDetrend)
    
    cfigopen('Navigator Phase Detrending')
    
    nframes = size(raw,1);
    ndat = size(raw,2);
    nleaves = size(raw,3);
    nslices = size(raw,4);
    ncoils = size(raw,5);
    
    plot(angle(raw(ceil(nframes/2),:,ceil(nleaves/2),ceil(nslices/2),ceil(ncoils/2))));
    hold on
    plot(PDfits(ceil(nframes/2),:,ceil(nleaves/2),ceil(nslices/2),ceil(ncoils/2)));
    plot(angle(raw_corr(ceil(nframes/2),:,ceil(nleaves/2),ceil(nslices/2),ceil(ncoils/2))));
    hold off
    
    legend('Uncorrected phase',sprintf('Polynomial fit (order %d)',orderPhaseDetrend),'Corrected phase',...
        'Location','southOutside');
    
    title(sprintf('Navigator phase detrending for example echo\n(frame: %d/%d, leaf: %d/%d, slice: %d/%d, coil: %d/%d)',...
        ceil(nframes/2),nframes,ceil(nleaves/2),nleaves,ceil(nslices/2),nslices,ceil(ncoils/2),ncoils));

end

function genfig_T2(raw,raw_corr,T2,T2curve)

    cfigopen('T2 filtering')
    
    nframes = size(raw,1);
    ndat = size(raw,2);
    nleaves = size(raw,3);
    nslices = size(raw,4);
    ncoils = size(raw,5);
    
    plot(reshape(abs(raw(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1))
    hold on
    plot(reshape(abs(raw_corr(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1))
    plot(reshape(abs(T2curve(ceil(nframes/2),:,ceil(nleaves/2),:,ceil(ncoils/2))),[],1))
    hold off
    
    legend('Uncorrected Echo','Corrected Echo','Fitted T2 Correction Curve')
    title(sprintf('T2 correction for frame %d, leaf %d, coil %d\nMean estimated T2: %.2fms', ...
        ceil(nframes/2),ceil(nleaves/2),ceil(ncoils/2),T2))
    
    drawnow
    
end

function genfig_psf(psf_3D,fov)

    cfigopen('Orthogonal Views of Point Spread Function')
    
    showov(log(abs(psf_3D)), 'Image (log)', 1, 1, ...
        fov(1)/2*[-1 1], ...
        fov(2)/2*[-1 1], ...
        fov(3)/2*[-1 1]);
    
    sgtitle('Point spread function')
    
    drawnow

end

function genfig_ovcoil(im,ksp,fov,Kmax,framen,nframes,coiln,ncoils)

    cfigopen('Orthogonal Views of Reconstruction for Current Coil')

    showov(abs(im), 'Image', 1, 2, ...
        fov(1)/2*[-1 1], ...
        fov(2)/2*[-1 1], ...
        fov(3)/2*[-1 1]);
    
    showov(log(abs(ksp)), 'K-space (log)', 2, 2, ...
        [sprintf("%.1f",-Kmax(1)) sprintf("%.1f",Kmax(1))], ...
        [sprintf("%.1f",-Kmax(2)) sprintf("%.1f",Kmax(2))], ...
        [sprintf("%.1f",-Kmax(3)) sprintf("%.1f",Kmax(3))]);

    sgtitle(sprintf('Reconstruction for frame %d/%d, coil %d/%d',...
        framen,nframes,coiln,ncoils))
    
    drawnow

end

function genfig_ovframe(im,ksp,fov,Kmax,framen,nframes)

    cfigopen('Orthogonal Views of Reconstruction for Current Frame')

    showov(abs(im), 'Image', 1, 2, ...
        fov(1)/2*[-1 1], ...
        fov(2)/2*[-1 1], ...
        fov(3)/2*[-1 1]);
    
    showov(log(abs(ksp)+1), 'K-space (log)', 2, 2, ...
        [sprintf("%.1f",-Kmax(1)) sprintf("%.1f",Kmax(1))], ...
        [sprintf("%.1f",-Kmax(2)) sprintf("%.1f",Kmax(2))], ...
        [sprintf("%.1f",-Kmax(3)) sprintf("%.1f",Kmax(3))]);

    sgtitle(sprintf('Reconstruction for frame %d/%d (all coils combined)',...
        framen,nframes))
    
    drawnow

end

function genfig_tSNR(tSNR,fov)
   
    cfigopen('tSNR Map')
    
    showov(abs(tSNR), 'Image', 1, 1, ...
        fov(1)/2*[-1 1], ...
        fov(2)/2*[-1 1], ...
        fov(3)/2*[-1 1]);
    
    sgtitle('Temporal Signal/Noise Ratio (tSNR)')
    
    drawnow

end
