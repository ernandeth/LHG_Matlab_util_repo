function [args, scaninfo, kinfo] = rt_rec_setup1(varargin)

if nargin == 0
    fprintf(stderr,'\n Usage: %s databuffer [OPTIONS] \n',mfilename);
    fprintf(stderr,'\n       OR         ');
    fprintf(stderr,'\n Usage: %s rawfile [OPTIONS]\n',mfilename);
    fprintf(stderr,'\n Recon Options\n');
    fprintf(stderr,'cp   Do conj phase recon\n');
    fprintf(stderr,'iter Do iterative (CG) recon\n');
    fprintf(stderr,'jnt  Do joint estimation recon\n');
    fprintf(stderr,'m    Generate field maps (for use with -h)\n');
    fprintf(stderr,'F #  Filter for smoothing fieldmap (default %d)\n',FSIZE);
    fprintf(stderr,'M #  Navigator correction for multishot (magnitude)\n');
    fprintf(stderr,'P #  Navigator correction for multishot (phase)\n');
    fprintf(stderr,'O #  Off-resonance demodulation (in Hz)\n');
    fprintf(stderr,'n #  Recon image size (default set by pulse sequence)\n');
    fprintf(stderr,'d #  Sample delay term (for calibration)\n');
    fprintf(stderr,'Q    Do out spiral data in in-out spiral acquisition\n');
    fprintf(stderr,'\nOutput Options\n');
    fprintf(stderr,'V    Read header only (no recon)\n');
    fprintf(stderr,'v    Verbose\n');
    fprintf(stderr,'l    Shift images to actual Rx-ed location\n');
    fprintf(stderr,'A    Output Analyze format (default)\n');
    fprintf(stderr,'sl   Output individual images\n');
    fprintf(stderr,'mat  Output matlab .mat files\n');
    fprintf(stderr,'com  Output complex data matlab .mat files\n');
    fprintf(stderr,'S    Output matlab k-space data (for diags)\n');
    fprintf(stderr,'a    Output matlab file for all coils in multicoil recon\n');
    fprintf(stderr,'R1   Rotation image 90 deg.\n');
    fprintf(stderr,'R2   Rotation image 180 deg. (also -R)\n');
    fprintf(stderr,'R3   Rotation image 270 deg.\n');
    fprintf(stderr,'fx   Flip images in x-direction\n');
    fprintf(stderr,'fy   Flip images in y-direction\n');
    fprintf(stderr,'q    Reverse slice ordering\n');
    fprintf(stderr,'x #  Shift images in x-direction by # pixels\n');
    fprintf(stderr,'y #  Shift images in y-direction by # pixels\n');
    fprintf(stderr,'z #  Zoom factor (default 1.0)\n');
    fprintf(stderr,'L    Do L2 norm with existing images (for spiral in/out)\n');
    fprintf(stderr,'t #  Recon temporal frame # only\n');
    fprintf(stderr,'s #  Recon slice # only (incompatible with -A)\n');
    fprintf(stderr,'o #  Offset output image numbers (adds # to temp frame)\n');
    fprintf(stderr,'I    For interleaved slice ordering\n');
    fprintf(stderr,'C #  To adjust output scaling (default = %f)\n',outscale);
    fprintf(stderr,'0    (zero) to output 0th time point (field map)\n');
    return
end

% set default values
args.pfile = 0;
args.rtype = 'ft'; % 'cp','ft','iter','jnt'
args.map = 0;  %  /* correction maps for hom. corr. */
args.numrot = 0; % /* extra rotations */
args.phs = 0; %     /* output phase images */
args.out = 0; % verbose
args.norec = 0;% , no recon
args.allcoils = 0; %/* output images for all coils */
args.L2 = 0; %      /* output L2 combintation of images */
args.imoffset = 0;
args.setnim = 0;
args.nim = 64;
args.iph = 0;
args.isl = 0;
args.samp1 = 0;
args.fsize = FSIZE;
args.lrshift = 0;
args.tbshift = 0;
args.zoomer = 1;
args.outscale = outscale;
args.pt = 0;
args.phasefact = 0;
args.svraw = 0; %   /* save mag raw data into images */
args.loc = 0; %     /* true locations, shift images */
args.rev = 0; %     /* reverse slice ordering */
args.outfmt = 'nifti'; %    /* NIFTI output format */
args.unf = 0; %     /* UNFOLD recon */
args.unfact = 1; %     /* UNFOLD recon */
args.zflg = 0; %    recon zero point
args.rsa = 0; %     /* Out flag in in-out spiral */
args.sliceorder = 1; %  /* 0=interleaved, 1=sequential */
args.flipx = 0;   % flip x axis at output
args.flipy = 0;   % flip x axis at output
args.complex = 0;


%
% read input args
%
argn = 1;
while (argn <= length(varargin));
    argtype = ( cell2mat(varargin(argn) ));
    argn = argn + 1;
    switch (argtype)
        case 'ft', args.rtype = 'ft';
        case 'cp', args.rtype = 'cp';
        case 'h', args.rtype = 'cp';
        case 'iter', args.rtype = 'it';
        case 'jnt', args.rtype = 'jn';
        case 'm',  args.map = 1;  %  /* correction maps for hom. corr. */
        case 'R1', args.numrot = 1; % /* extra rotations */
        case 'R2', args.numrot = 2; % /* extra rotations */
        case 'R',  args.numrot = 2; % /* extra rotations */
        case 'R3', args.numrot = 3; % /* extra rotations */
        case 'v', args.out = 1; % verbose
        case 'V', args.out = 1; args.norec = 1;% verbose, no recon
        case 'a', args.allcoils = 1; %/* output images for all coils */
        case 'L', args.L2 = 1; %      /* output L2 combintation of images */
        case 'o', args.imoffset = getnum( cell2mat(varargin(argn) )); argn = argn+1;
            %		 /* offset in output image numbers */
        case 'n', args.setnim = 1;
            args.nim = getnum(cell2mat(varargin(argn)));
            argn = argn+1;  %  image size
        case 't', args.iph = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 's', args.isl = getnum(cell2mat(varargin(argn)));  argn = argn+1;
        case 'd', args.samp1 = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'F', args.fsize = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'x', args.lrshift = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'y', args.tbshift = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'z', args.zoomer = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'C', args.outscale = getnum(cell2mat(varargin(argn))); argn = argn+1;
            %      case 'O', args.pt = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'P', args.phasefact = getnum(cell2mat(varargin(argn))); argn = argn+1;
        case 'S', args.svraw = 1; %   /* save mag raw data into images */
        case 'l', args.loc = 1; %     /* ture locations, shift images */
        case 'q', args.rev = 1; %     /* reverse slice ordering */
        case 'A', args.outfmt = 'anal'; %    /* ANALYZE output format */
        case 'N', args.outfmt = 'nifti'; %    /* ANALYZE output format */
        case 'sl', args.outfmt = 'sl'; %
        case 'mat', args.outfmt = 'mat'; %
        case 'com', args.complex = 1;
        case 'U', args.unf = 1; %     /* UNFOLD recon */
        case 'Q', args.rsa = 1; %     /* Out flag in in-out spiral */
        case '0', args.zflg = 1; %    recon zero time point
        case 'I', args.sliceorder = 0; %  /* 0=interleaved, 1=sequential */
        case 'fx', args.flipx = 1;   % flip x axis at output
        case 'fy', args.flipy = 1;   % flip x axis at output
        otherwise,
            fprintf(stderr,'%s: Did not recognize option %s\n',mfilename,argtype);
    end % switch
end % while
