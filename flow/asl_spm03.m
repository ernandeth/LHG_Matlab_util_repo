function result = asl_spm03(args)
%
% function result = asl_spm02(args)
%
% (c) 2010-2019 Luis Hernandez-Garcia @ University of Michigan
% report problems to hernan@umich.edu
%
% This function processes ASL data from a raw k-pace file to a
% QUANTITATIVE PERFUSION ACTIVATION map of the deisred effects
%
% some notes:  the name of the file that the work is being done is stored
% in the variable "workFile"
%
% if you call the function without arguments,
% you get the args structure with the default values:
% and the program exists.
%
%
%     % default values:

if nargin<1
    % if you call the function without arguments,
    % you get the args structure with the default values:
    % and the program exists.
    %
    % default values:
    args.inFile = [];
    args.doDespike = 0;
    args.doRecon=0;
    args.doZgrappa= 0;
    args.doSliceTime=0;
    
    args.doRealign = 0;
    args.smoothSize= 8;
    args.subType = 0;
    args.physCorr = 0;
    args.physFile= [];
    args.CompCorr = 0;
    
    args.anat_img = [];
    args.template_img = [];
    args.spat_norm_series = 1;
    
    args.doQuant = 0;
    args.Diff_img = 'mean_sub.img';  
    args.SpinDens_img = 'SpinDensity.img';
 
    args.aslParms.TR = 3;
    args.aslParms.Ttag = 1.4;
    args.aslParms.Tdelay = 1.2;
    args.aslParms.Ttransit = 1.2;
    args.aslParms.inv_alpha = 0.8;
    args.aslParms.disdaqs = 2;
    
    
    args.doGLM = 0;
    args.designMat = [];
    args.isSubtracted = 1;
    args.contrasts = [0 0 -1 1];
    args.doQuant_GLM = 0;
    args.BaseFlow_img = 'wFlow.img';
    args.doGlobalMean = 0;
    args.doLightbox = 0;
    args.doOrtho = 1;
    args.FlowScaleFactor = 1;
    args.subOrder = 1;
    
    result = args;
    return
end
%
% save arguments for future use
save asl_spm03_params.mat args

try
    
%% First figure out the input working file name

[pth name ext] = fileparts(args.inFile(1,:));

if strcmp(ext, '.nii')
    mkdir('Merged')
    cd('Merged')
    
    fprintf('\nNow working in directory: \n%s\n',pwd)
    
    disp('Now Merging Files')
    
    spm_file_merge(args.inFile, [pwd '/merged.nii']);
    
    workFile = [pwd '/merged.nii'];

end

workFile = args.inFile(1,:);
[pth name ext] = fileparts(workFile);
if isempty(pth)
    pth=pwd;
end
cd(pth)
fprintf('\nNow working in directory: \n%s\n',pth);
if strcmp(ext, '.img');
    fprintf('\nConverting %s from Analyze to NIFTI ...\n', workFile);
    [d h]=read_img(workFile);
    h = avw2nii_hdr(h);
    write_nii( [name '.nii'], d, h, 0);
    workFile = fullfile(pth, [name '.nii'])
    volFile = workFile;
end

%% recon section .. if found white pixel arttefact -> use despiker
if args.doRecon
    
    if args.doDespike
        fprintf('\ndespiking .... %s\n', workFile);
        despiker_ASL(workFile, 4, 0);
        workFile = ['f_' workFile];
    end
    
    %{
    if args.doRecon==1
        fprintf('\ndoing 2D recon on ....%s\n', workFile);
        sprec1(workFile, 'l', 'fy','N');
    end
    %}
    fprintf('\ndoing 3D recon on ....%s\n', workFile);
    fprintf('\nwarning ... REMOVING old *VOL* FILES FROM DIRECTORY '); pause(3)
    !rm *vol*
    
    % include Z grappa option for recon
    if args.doZgrappa == 1,
        fprintf('\n(Using 1-D GRAPPA along Z axis)');
        sprec1_3d_grappaz(workFile, 'l', 'fy','N', 'C', 1, 'grappaz', args.M0frames);
        
        fprintf('\nSplitting the time series into:');
        fprintf('\n    - mean of %d Calibration Frames (spin density)', args.M0frames);
        fprintf('\n    - BGS suppressed label-control frames ...'); 
        load fullysampled.mat
        tmp = dir('vol*.nii');
        h = read_nii_hdr(tmp(1).name);
        h = nii2avw_hdr(h);
        h.tdim = 1;
        write_img('SpinDensity.img', fullimg(:),h);
    else
        fprintf('\n(No GRAPPA along Z axis)');
        sprec1_3d_grappaz(workFile, 'l', 'fy','N', 'C', 1, 'l');
        
        if args.M0frames>0
            fprintf('\nSplitting the time series into:');
            fprintf('\n    - mean of %d Calibration Frames (spin density)', args.M0frames);
            fprintf('\n    - BGS suppressed label-control frames ...');
            
            tmp = dir('vol*.nii');
            workFile = tmp(1).name;
            volFile = workFile;
            
            [dat h2] = read_nii_img(volFile);
            m0 = dat(1:args.M0frames, :);
            m0 = mean(m0,1); 
            dat = dat(args.M0frames+1:end, :);
            
            h2.dim(5) = h2.dim(5) - args.M0frames;
            write_nii(volFile, dat,h2, 0);
            
            h = nii2avw_hdr(h2);
            h.tdim = 1;
            write_img('SpinDensity.img', m0(:),h);
        end
        
    end
    
    tmp = dir('vol*.nii');
    workFile = tmp(1).name;
    volFile = workFile;
    
    
end

%% Slice Timing correction for ASL
if args.doSliceTime
    
    fprintf('\ndoing slice timing on ....%s\n', workFile);
    asl_sliceTimer(workFile, args.aslParms.TR, args.aslParms.Ttag+args.aslParms.Tdelay);
    workFile = ['a' workFile];
    
end



%% REALIGN images
if args.doRealign
    %  realignment with MCFLIRT
    %  !setenv FSLOUTPUTTYPE NIFTI
    %  str = sprintf('!mcflirt -in %s -out rvol -refvol 0 -cost normcorr -verbose 1 -stats -plots -mats', n);
    %  eval(str)
    fprintf('\ndoing realignment on ....%s\n', workFile);
    % realignment with SPM
    opts.rtm=1;
    spm_realign(workFile, opts);
    spm_reslice(workFile);
      
    [pth, root, ext] = fileparts(workFile);
    workFile = fullfile(pth, ['r' root ext]);
    
end

%%  Gaussian Smoothing
if args.smoothSize >0
    
    
    sz=args.smoothSize;
    
    fprintf('\nsmoothing  ....%s with %f kernel \n', workFile, sz);

    % spm_smooth doesn't handle paths gracefully:
    [pth, root, ext] = fileparts(workFile);
    curDir = pwd;
    if ~isempty(pth); cd(pth); end
    workFile = [root ext];
    spm_smooth(workFile,['s' workFile],[ sz sz sz], 4 );
    workFile = ['s' workFile];
    cd(curDir)
end

%%  Physio Correction
if args.physCorr==1
    
    fprintf('\ndoing RETROICOR on ....%s (must be unsubtracted)\n', workFile);
    if args.doSliceTime==1 , timeType=1, else timeType=0; end
    
    [d h] = read_nii_img(workFile);
    Nslices = h.dim(4);
    args.aslParms.disdaq = 2;
    
    % step 1: read the physio  data from the scanner and change its format
    % this generates a file called physio.dat
    physdata = convertEXphysio(args.physFile, 0.025);
    
    % step 2: create a matrix with the physio data made up of
    % basis functions
    PhysioMat = mkASLPhysioMat(...
        'physio.dat',...
        0.025, ...
        args.aslParms.disdaq, ...
        Nslices, ...
        args.aslParms.TR, ...
        args.aslParms.Ttag + args.aslParms.Tdelay,...
        timeType);
    
    % step 3:  estimate the parameters of the design matrix and remove iit from
    % the data .  This will generate a 4D image file called residuals.nii.
    rmReg(workFile, PhysioMat, 2);
    
    workFile = ['residuals.nii'];
end

%%
if args.subType > 0
    ortho_args = ortho2005;
    ortho_args.tseries_file =  workFile;
    ortho_args.ROItype = 'sphere';
    ortho_args.ROIsize = 2;
    ortho_args.doMovie = 0;
    ortho_args.interact = 0;
    
    ortho2005(ortho_args);
    title(sprintf('UN - Subtracted Time Series'))
    set(gcf,'Name', 'OrthoView of Raw Images')
    close(findobj('Name','Histograms'))
    
    % subtraction section
    h = read_nii_hdr(workFile);
    Nframes = h.dim(5)
    warning off
    tfile = workFile;
    
    switch args.subType
        
        case 1
            fprintf('\ndoing pairwise subtraction on ...%s\n', workFile);
            !rm sub.img sub.hdr
            [p, rootname,e] = fileparts(workFile);
            aslsub(rootname, 1, 1, Nframes, 0, args.subOrder, 0);
            
            
            ms = read_img('mean_sub');
            if sum(ms(:)) < 0
                fprintf('\n WARNING:  reversing the subtraction order! \n')
                aslsub(rootname, 1,  1, Nframes, 0, ~(args.subOrder), 0);
                ms = lightbox('mean_sub',[-200 200],[]);
            end
            workFile = 'sub.img';
            title('Mean subtraction (pairwise)');
            set(gcf,'Name', 'Slice View of Mean Difference Image')
            
        case 2
            fprintf('\ndoing surround subtraction on ...%s\n', workFile);
            !rm sub.img sub.hdr
            [p, rootname,e] = fileparts(workFile);
            aslsub_sur(rootname,  1, Nframes, 0, args.subOrder);
            
            ms = read_img('mean_sub');
            if sum(ms(:)) < 0
                fprintf('\n WARNING:  reversing the subtraction order! \n')
                aslsub_sur(rootname, 1,Nframes, 0, ~(args.subOrder));
            end
            workFile = 'sub.img';

%             figure
%             ms = lightbox('mean_sub',[],[]);
%             p = get(gcf, 'Position');
%             set(gcf,'Position', p + [3 -2 0 0]*100);
%             title('Mean Subtraction (surround)');
%             set(gcf,'Name', 'Slice View of Mean Difference Image')
%             colormap hot
    end
    %
    %

  

    ortho_args = ortho2005;
    ortho_args.anat_file =  'mean_sub';
    ortho_args.tseries_file =  workFile;
    ortho_args.wscale = [0 200];
    ortho_args.ROItype = 'sphere';
    ortho_args.ROIsize = 2;
    ortho_args.doMovie = 0;
    ortho_args.interact = 0;
    ortho2005(ortho_args);
      

    colormap hot
    p = get(gcf, 'Position');
    set(gcf,'Position', p + [2 -1 0 0]*100);
    title(sprintf('Subtracted Time Series'))
    set(gcf,'Name', 'OrthoView of Mean Difference Image')
    

    close(findobj('Name','Histograms')) 

end
%% Physio correction section using CompCor:
%  use this with subtracted data only
if args.CompCorr==1
    fprintf('\nPerforming PCA CompCorr on %s  ...\n', workFile)
    [dirty hdr] = read_img(workFile);
    
    
    [clean junkcoms] = compcor12(dirty, hdr, 12);
            
    fprintf('\nWriting ....clean_%s version \n', workFile);
    [pth root ext] = fileparts(workFile);
    tmp = ['clean_' root ext];
    
    if isfield(hdr, 'tdim')
        write_img(tmp, clean, hdr);
    else
        write_nii(tmp, clean, hdr, 0);
    end
    
    if ~isempty(args.designMat)
        % decorrelate the designmatrix out of the confounds
        X = args.designMat; % regressors of interest
        pr = pinv(X);
        fprintf('\nDecorrelating junk regressors from CompCor  ...\n')
        whos junkcoms X
        for n=1:size(junkcoms,2)
            reg = junkcoms(:,n);
            reg = reg - X*(pr*reg);
            % mean center the confounds
            reg = reg -mean(reg);
            junkcoms(:,n) = reg;
        end
        
        args.designMat = [args.designMat junkcoms];
        args.contrasts = [args.contrasts zeros(size(args.contrasts,1),size(junkcoms,2))] ;
        x = args.designMat ;
        

        fprintf('\nWriting out adjusted (compcor) design matrix : designMat_compcor.dat ...\n')
        save designMat_compcor.dat x -ascii
        
        figure; imagesc(args.designMat);            
        set(gcf,'Name', 'Compcor Identified noise patterns in the time series ')
        title('Design Matrix PLUS CompCor confounds');drawnow
    end

end
%%
if args.doQuant==1
    
    fprintf('\nCalculating mean baseline CBF from consensus paper model ...');
    
    M0frames = args.M0frames;  % the first frames do not have background suppression
    inv_alpha = args.inv_alpha;
    flip = 30 * pi/180;
    Ttag = args.Ttag;
    TR = args.TR;
    pid = args.Tdelay;
    T1 = args.T1;
    
    % Check to see if this is the vasc file from the GE scanner
    % they have two time frames: difference and spin density images.
    % If that's the case, we'll split it into the 
    % INDIVIDUAL subtraction and spin density image files 
    % otherwise, we ignore the inFile for the CBF calculation
    
    %Vo = spm_file_split(args.inFile);
    %if length(Vo)==2
    %    movefile(Vo(2).fname, 'SpinDensity.nii', 'f');
    %    movefile(Vo(1).fname, 'mean_sub.nii', 'f');
    %    
    %    args.Diff_img = 'mean_sub.nii';
    %    args.SpinDens_img = 'SpinDensity.nii';
    %end
    
    [msk h] = read_img(args.SpinDens_img);
    thresh = 1*std(msk(:));
    msk(msk<=thresh) = 0;
    msk(msk>0) = 1;
    
    f = calc_cbf_wp(args.Diff_img, args.SpinDens_img, inv_alpha, Ttag, TR, pid, T1);
    
    [d h] = read_img('Flow.img');
    d = msk .* d * args.FlowScaleFactor;
    write_img('Flow.img',d,h);
    
    % Read Flow.img and multiply by scale factor here
    
    figure
    subplot(221),  lightbox('mean_sub');  title('mean subtraction')
    subplot(222),  lightbox('sSpinDensity'); title('smoothed spin density')
    subplot(223),  lightbox('Flow', [0 60], []) ; title('Perfusion (ml/min/100g)')
    
end

%% coregistration section:
if ~isempty(args.anat_img)
    
    fprintf('\nCoregistering Structural to mean subtraction image ...');
    flags.cost_fun='nmi';
    flags.tol = [0.01 0.01 0.01 0.001 0.001 0.001];
    
    Vref = spm_vol(args.Diff_img);
    Vtgt = spm_vol(args.anat_img);
    
    x = spm_coreg(Vref,Vtgt, flags);
    
    % set the affine xformation for the output image's header
    mat = spm_matrix(x);
    xform = inv(mat)*Vtgt.mat ;
    
    spm_get_space(Vtgt.fname,xform);
    
    
    % normalization fMRI to template
    if ~isempty(args.template_img)
        fprintf('\nSpatially Normalising Structural to template image ...');
        
        Vref = spm_vol(args.template_img);
        Vtgt = spm_vol(args.anat_img);
        
        spm_normalise(Vref, Vtgt, 'mynorm_parms.mat');
        
        if args.spat_norm_series==1
            % apply the normalization to the sub images
            fprintf('\nApplying Normalization to time series %s ...', workFile);
            h = read_nii_hdr('sub.hdr');
            [pth root ext] = fileparts(workFile)
            workFile = fullfile(pth , [root '.img']);
            
            for n=1:h.dim(5)
                subfiles{n} = [workFile ',' num2str(n)];
                spm_write_sn( subfiles{n}  , 'mynorm_parms.mat');
            end
            workFile = fullfile(pth , ['w' root '.img']);
        end
        
        fprintf('\nApplying Normalization to anatomical and mean_sub ...');
        spm_write_sn( './mean_sub.img'  , 'mynorm_parms.mat');
        spm_write_sn( args.anat_img , 'mynorm_parms.mat');
        
        if args.doQuant
            fprintf('\nApplying Normalization to quantification images ...');
            spm_write_sn( './Flow.img'  , 'mynorm_parms.mat');
            spm_write_sn( './sSpinDensity.img'  , 'mynorm_parms.mat');
        end
        
    end
end

%% calculate a mean global at every image to create a global mean confound
if args.doGlobalMean==1
    fprintf('\nCalculating the Global mean at every time point ....%s\n', workFile);
    [tmp h] = read_img(workFile);
    gm = mean(tmp,2);
    gm = gm - mean(gm);
    gm = gm / max(gm);
    
    fprintf('\nIncluding the Global mean into the Design Matrix %s ....\n', workFile);
    args.designMat = [args.designMat gm];
    args.contrasts = [args.contrasts zeros(size(args.contrasts,1),1) ] ;
end

%% Rescaling time series so that the spatial-temporal mean is 1000 (inside mask);
if args.doGlobalMean==2
    
    fprintf('\nRescaling time series so that the spatial-temporal mean is 1000 (inside mask) ...');
    % the mask is set to all the pixelswhose signal is above one standard deviation
    global_scale(workFile,1);
    
end

%% Estimation of Parameters and Statistical Maps
if args.doGLM
    
    fprintf('\nEstimating GLM on %s ....\n', workFile);
    
    figure
    subplot(121)
    imagesc(args.designMat);
    title('This is the design matrix')
    subplot(222)
    imagesc(args.contrasts);
    title('These are the contrasts')
    colormap gray
    
    if args.isSubtracted==0
        flags.doWhiten=1;
    else
        flags.doWhiten=0;
    end
    
    flags.doWhiten=0;
    flags.header=[];
    
    spmJr(workFile, args.designMat ,args.contrasts, flags);
    
    % Write out the thresholded percent signal changes
    contrast2percent(1000, 2.5);
    
    % show the different Zmaps as a diagnostic:
%     for z=1:size(args.contrasts,1);
%         figure; lightbox(sprintf('Zmap_%04d.img', z),  [], 4);
%     end
end



%% convert betas to perfusions!
if args.doQuant_GLM
    % If the template was specified, this means that we did
    % a normalization on the images.

    bf_name = args.BaseFlow_img;
    fprintf('\n\nQuantifying flow contrasts.  Baseline is:  %s', bf_name);  
    [bf hf] = read_img(bf_name);
    pct_images = dir('percent*.img');
    
    for p=1:length(pct_images)
        [pct h] = read_img(pct_images(p).name);
        conFlow = pct .* bf / 1000 / 100;  % the percent images were scaled !
        conFlow(isnan(conFlow)) = 0;        
        write_img(sprintf('ConFlow_%04d.img',p), conFlow, h);
    end
        
end

%% making nice overlays of the flows
% if args.doLightbox ==1
%     
%     [flows h] = read_img('ExpFlows');
%     flows = reshape(flows, h.tdim,  h.xdim*h.ydim*h.zdim);
%     
%     f0 = reshape(flows(1,:), h.xdim,h.ydim,h.zdim);
%     
%     for f=1:size(flows,1)
%         f_act = reshape(flows(f , :), h.xdim,h.ydim,h.zdim) ;
%         figure
%         act_lightbox( f0, f_act, [10 80], [1 25], 5);
%         title(sprintf('Contrast number %d', f));
%     end
% end

%% DIsplay the activation maps
if args.doLightbox ==1 && args.doQuant_GLM
    bf_name = args.BaseFlow_img;
    if args.spat_norm_series==1
       bf_name = 'wmean_sub.img';
    end
    [underlay h] = read_img(bf_name);
    underlay = reshape(underlay,[h.xdim h.ydim h.zdim]);

    
    for f=1:size(args.contrasts,1)
        
        [zmap h] = read_img(sprintf('conFlow_%04d',f));
        zmap = reshape(zmap,[h.xdim h.ydim h.zdim]);
        figure
        act_lightbox( underlay, zmap, [0 2e2], [5 50], [6]);
        title(sprintf('Perfusion Contrast number %d: [%s]', f, num2str(args.contrasts(f,:))));
        set(gcf,'Name','Perfusion Change overlaid on Baseline Perfusion');
        ylabel('Perfusion')
        p = get(gcf, 'Position');
        set(gcf,'Position', p+ [50+(2+f) (2*f) 1 1]*100);
    end
end

%% DIsplay the activation maps
if args.doLightbox ==1
    bf_name = args.BaseFlow_img;
    if args.spat_norm_series==1
       bf_name = 'wmean_sub.img';
    end
    [underlay h] = read_img(bf_name);
    underlay = reshape(underlay,[h.xdim h.ydim h.zdim]);

    
    for f=1:4 %size(args.contrasts,1)
        
        [zmap h] = read_img(sprintf('Zmap_%04d',f));
        zmap = reshape(zmap,[h.xdim h.ydim h.zdim]);
        figure
        act_lightbox( underlay, zmap, [0 2e2], [2.5 10], [6]);
        title(sprintf('Contrast number %d: [%s]', f, num2str(args.contrasts(f,:))));
        set(gcf,'Name','Activation map (Z) overlaid on difference image');
        ylabel('Z score')
        p = get(gcf, 'Position');
        set(gcf,'Position', p+ [(2+f) (2*f) 1 1]*200);

    end
end

%% if needed, we can put it in ortho to explore the time course
if args.doOrtho
    
    tfile=workFile;
    contrasts = args.contrasts
    if args.subType
        tfile='sub.img';
    end
    fprintf('\nDisplaying the last contrast in Orthogonal views');
    
    asl_args = args;
    clear global args
    
    ortho2005([],...
        'anat_file', 'Bhats.img', ...
        'tseries_file', tfile, ...
        'spm_file', sprintf('Zmap_%04d.img', size(contrasts,1)),...
        'ROItype', 'sphere', ...
        'ROIsize', 2, ...
        'threshold', 2 ...
        );
   close(findobj('Name','Histograms')) 
end
%
%% Adding new stuff below (8/21/13)
%
% Make a mean image of the functional(s)
% [data h] = read_img(workFile);
% [p n e]=fileparts(workFile);
% if e ~='.nii'
%     h=avw2nii_hdr(h);
% end
% meandata = mean(data, 1);
% h.dim(5)=1;
%write_nii('mean_func.nii', meandata, h,0);


%% Preprocessing of anatomical files- do this if files are not already reconstructed.

% if args.doBuildAnatomy
%     curDir=pwd;
%     cd (args.anatomyDir);
%     % !buildanatomy;   % not working right now
%     str =['!cp eht1* ' curDir];
%     eval(str)
%     str =['!cp ht1* ' curDir];
%     eval(str)
%     str =['!cp t1* ' curDir];
%     eval(str)
%
%     cd(curDir)
% end

%% Coregistration of functionals and statistical maps.
% if args.doCoreg % && args.useSPGR
%     x=spm_coreg(args.overlayfile, 'mean_sub.img' );
%     M=(spm_matrix(x));
%
%     tmpM = spm_get_space(args.overlayfile);
%     spm_get_space(args.overlayfile, M*tmpM);
%     spm_get_space(args.overlayfile(3:end), M*tmpM);
%
%
%     x=spm_coreg(args.spgrfile,args.overlayfile);
%     M=(spm_matrix(x));
%     tmpM=spm_get_space(args.spgrfile);
%     spm_get_space(args.spgrfile, M*tmpM);
%     spm_get_space(args.spgrfile(3:end), M*tmpM);
%
%
%     name_cells = {...
%         args.spgrfile, ...
%         args.overlayfile,...
%         args.spgrfile(3:end), ...
%         args.overlayfile(3:end),...
%         'mean_sub.img', ...
%         'sub.img'};
%
%     zmaps=dir('Zmap_*.img');
%     for n=1:length(zmaps)
%         name_cells{n+4} = zmaps(n).name;
%     end
%
%     %spm_reslice(name_cells, struct('which',1,'mask',0,'mean',0));
%
%     zmaps=dir('ConBhat_*.img');
%     for m=1:length(zmaps)
%         name_cells{m+n+4} = zmaps(m).name;
%     end
%
%     spm_reslice(name_cells, struct('which',1,'mask',0,'mean',0));
% end
%
% %% DO the spatial normalization
% if isfield(args, 'norm_ref')
%     if ~isempty(args.norm_ref)
%
%         spm_normalise('/export/home/hernan/matlab/spm8/templates/T1.nii',...
%             args.norm_ref, ...
%             'mynorm_parms');
%
%         for n=1:length(args.norm_list)
%             spm_write_sn(args.norm_list{n}, 'mynorm_parms.mat');
%         end
%
%     end
% end

catch asl_spm_errors
    asl_spm_errors
    fprintf('\n\nErrors happened ... Exiting. ');
    fprintf('\nThe error structure is in "asl_spm_errors.mat":');
    save asl_spm_errors.mat asl_spm_errors
    return
end

return

