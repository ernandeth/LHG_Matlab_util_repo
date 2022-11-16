 function sprec2(pfile,varargin)
%function sprec2(pfile,varargin)
% Copyright 2005
% Douglas C. Noll and others 
% University of Michigan 
%
% call sprec2 with no arguments to find options
% 
% to allow immeidate access to matlab recon variables
% edit this file to comment out the top line and change
% scriptmode = 1;
%

scriptmode = false;     % set to 'true' if running sprec2 in matlab
global REC_PATH;        % path to directory containing recon

% setup path
% NOTE: edit path to point to the
% directory where sprec2.m is stored
if scriptmode && isempty(REC_PATH)
    REC_PATH = '~fmrilab/recon/matlab/v2';
end
setup_path(REC_PATH);

if scriptmode
  % in script mode, this is how recon is setup
    pfile_dir = strcat(REC_PATH, '/test');
    pfile = strcat(pfile_dir,'/070915_results/070915_P37888.7')
%    pfile = strcat(pfile_dir,'/070801ss_results/070801_P08192.7')
%    pfile = strcat(pfile_dir,'P44032.7')
%    pfile = strcat(pfile_dir,'P47104.7')
%    pfile = strcat(pfile_dir,'P06656.7')
%    pfile = strcat(pfile_dir,'P12800.7')
%    pfile = strcat(pfile_dir,'P16384.7')
%    pfile = strcat(pfile_dir,'P20480_splx3rf_10temp.7')
%    pfile = strcat(pfile_dir,'P22528_splx3rf_10temp_in_out.7')
%    pfile = strcat(pfile_dir,'P46080.7')
%    pfile = strcat(pfile_dir,'P52224.7')

% test files from CY
    pfile = '~chunyuy/MRdata/10Aug07/P29184.7'  % 4-shot
%    pfile = '~chunyuy/MRdata/10Aug07/P30208.7'  % 16-shot, 256x256
%    pfile = '~chunyuy/MRdata/10Aug07/P30720.7'

% other random test files
%    pfile = '/net/burton/data/volafsso/Data/BOLDapproach/041213ro_rw/ragnar/Pragnar_approach'
%    pfile = '/net/ed/data/grlee/flow/031211/P02560.7'
%    pfile = '/net/ed/home/volafsso/research/dR2/src/real_data/P15360.7'
%    pfile = '/net/burton/data/volafsso/Data/daqdel/041022recal/Psplx2_ddo-11'

% read p-file header
    [args,scaninfo,kinfo] = rec_setup1(pfile,'it','h','v','n',64,'t',1,'l','fy','fx');%,'t',6);
% uncomment to estimate field map
%    [args,scaninfo,kinfo] = rec_setup1(pfile,'it','m','v','n',64,'l','fy','fx');%,'t',6);
else
    if length(pfile) == 0
        rec_setup1;
        return
    else
        [args,scaninfo,kinfo] = rec_setup1(pfile,varargin{:});
    end
end % if scriptmode == 1

if args.out == 1
    args
    scaninfo
end
if args.norec == 1
    return
end

if str2num(scaninfo.scan_date(end-2:end)) <= 104
    fid = fopen(args.pfile,'r','b');
else
    fid = fopen(args.pfile,'r','l');
end

nim = args.nim;
iph = args.iph;
isl = args.isl;
zflg = args.zflg;
mapdda = scaninfo.mapdda;
unfact = args.unfact;
nphases = scaninfo.nphases;
nslices = scaninfo.nslices;
ncoils = scaninfo.ncoils;
npr = scaninfo.npr;
concat = scaninfo.concat;
datestr = scaninfo.scan_date(1:9);
datestr(3) = '_';
datestr(6) = '_';

if strcmp(args.outfmt,'nifti')
    pixperframe=nim*nim*nslices;
    nii_start = 1;
    nii_stop = pixperframe;
    niidata = zeros(pixperframe * (nphases-1),1);
end


if (args.rtype == 'cp') | (args.rtype == 'ft')
    % for gridding recon
    % set up kaiser-bessel functions
    scfactor = 0.04;  % arbitrary output scaning factor (depends on kb kernel, etc.)
    OVER = 2; % oversampling factor
    NW = 1024; % number of samples of k-b kernel
    kb_gridl = 1.5; % gridding width (in orig coords)
    kb_beta = 6.7; % Jackson's betas 1.5 = 6.7, 2.0 = 9.1, 2.5 = 11.5
    kbin = [-NW/2:NW/2]/NW;
    kbweights = besseli(0,kb_beta*sqrt(1.0-(2.0*kbin).^2));
    kx = []; ky = []; tt = [];
    for prnum = 1:npr
        kx = [kx real(kinfo.k.*exp(-i*2*pi*(prnum-1)/npr))*OVER];
        ky = [ky imag(kinfo.k.*exp(-i*2*pi*(prnum-1)/npr))*OVER];
        tt = [tt kinfo.t];
    end % for prnum = 1:npr
    ndat = length(kinfo.t);
    %r /* kaiser-bessel correction */
    wcmax = sinh(sqrt(kb_beta^2))/sqrt(kb_beta^2);
    wctmp = (2*pi*kb_gridl*[-nim/2:nim/2-1]/nim/OVER).^2 - kb_beta^2 + eps;
    wc1d = real(sqrt(wctmp)./sin(sqrt(wctmp))*wcmax);
    wc = wc1d'*wc1d;
    [xx yy] = meshgrid([-nim/2:nim/2-1]);
    wc = wc.*(abs(xx + i.*yy) < 0.55*nim)*scfactor*args.outscale;
    subim = [-nim/2:nim/2-1] + nim*OVER/2 + 2;
    % set up apodization
    apod = fermi(nim*OVER,max(abs(kinfo.k))*OVER,-1);
    if ((args.rtype == 'cp')|(args.map == 1))
        % set up standard temporal interpolators
        tmax = max(kinfo.t);
        % for now, make the max segment length 4 ms long - make this user selectable
        L = ceil(tmax/0.004);
        tau = tmax/L;
        t2 = [0:L]*tau;
        % for now, make the temporal interpolators min-max interpolators assuming
        % a generic triangular histogram with a FWHM of 100 Hz
        df = 100;
        [ll llt] = meshgrid([0:L]);
        ggt2 = sinc((ll-llt)*tau*df).^2;
        gb2 = sinc((([0:L]'*tau*ones(size(kinfo.t)) - ones([L+1 1])*kinfo.t)*df)).^2;
        beta = 0; %.00001;
        % a2 has the 1D temporal interpolators
        a2 = inv(ggt2 + beta*eye(L+1))*gb2;
        % tints has the gridded versions of a2
        % make time map
        datgr = reshape(kbgrid_mex(kx,ky,tt,ones(size(kx)),nim*OVER,kb_gridl*OVER,kbweights(:)),[nim nim]*OVER);
        tmap = real(datgr)./(imag(datgr)+eps);
        tints = zeros(nim*OVER,nim*OVER,L+1);
        for l = 0:L
            tmpint = interp1(kinfo.t,a2(l+1,:),tmap);
            tints(:,:,l+1) = tmpint.*apod;
        end
        tints(find(isnan(tints))) = 0;
    end
    kdens = kinfo.kdens.*kinfo.kmod;
    [kkx kky] = meshgrid([0:nim*OVER-1]);
    fftshmask = (-1).^(kkx + kky);
end

if strcmp(args.rtype,'it') % for iter
    kx = []; ky = []; tt = [];
    for prnum = 1:npr
        kx = [kx real(kinfo.k.*exp(-i*2*pi*(prnum-1)/npr))];
        ky = [ky imag(kinfo.k.*exp(-i*2*pi*(prnum-1)/npr))];
        tt = [tt kinfo.t];
    end % for prnum = 1:npr
    ndat = length(kinfo.t);
    [xx yy] = meshgrid([-nim/2:nim/2-1]);
    immask = (abs(xx + i.*yy) <= nim/2);
end


if (args.map == 1) % do field maps
    if (args.rtype == 'cp')
        rtype = 'ft';
        fm = zeros(nim,nim,nslices);
        stsl = 0; ensl = nslices; 
        stph = 0; enph = 2; iph = 0; args.unfact = 1; args.unf = 0;
        for slnum = stsl:ensl-1
            for phnum = 0:1
                % do field map recon here
                dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
                immap = zeros([nim nim]);
                for coilnum = 1:ncoils 
                    % do standard recon
                    docprec;
                      
                    if (phnum == 0)
                        imkbsv = imkb;
                        if (scaninfo.revflg == 2)
                            if (im1nz == 1)
                                im1sv = im1;
                            end
                        end
                    else 
                        imkb = imkb .* conj(imkbsv);
                        if (scaninfo.revflg == 2)
                            if (im1nz == 1)
                                imkb = imkb + im1 .* conj(im1sv);
                            end
                        end
                        immap = immap + imkb;
                    end % if (phnum == 0)
                end % for coilnum
            end % for phnum
            % add fancy smoother
            mnfm = mean(immap(:));
            fmtmp = angle(immap.*exp(-i*angle(mnfm)));
            fmmag = sqrt(abs(immap));
            fmsm = angle(conv2(fmmag.*exp(i*fmtmp),ones(args.fsize),'same'))+angle(mnfm);
            fm(:,:,slnum+1) = fmsm/(scaninfo.mapdel*1e-6)/2/pi; % in Hz
        end % for slnum

        % do it again using cp
        rtype = 'cp';
        for slnum = stsl:ensl-1
            for phnum = 0:1
                % do field map recon here
                dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
                immap = zeros([nim nim]);
                for coilnum = 1:ncoils 
                    % do standard recon
                    docprec;
      
                    if (phnum == 0)
                        imkbsv = imkb;
                        if (scaninfo.revflg == 2)
                            if (im1nz == 1)
                                im1sv = im1;
                            end
                        end
                    else 
                        imkb = imkb .* conj(imkbsv);
                        if (scaninfo.revflg == 2)
                            if (im1nz == 1)
                                imkb = imkb + im1 .* conj(im1sv);
                            end
                        end
                        immap = immap + imkb;
                    end % if (phnum == 0)
                end % for coilnum
            end % for phnum
            % add fancy smoother
            mnfm = mean(immap(:));
            fmtmp = angle(immap.*exp(-i*angle(mnfm)));
            fmmag = sqrt(abs(immap));
            fmsm = angle(conv2(fmmag.*exp(i*fmtmp),ones(args.fsize),'same'))+angle(mnfm);
            fm(:,:,slnum+1) = fmsm/(scaninfo.mapdel*1e-6)/2/pi; % in Hz
        end % for slnum
    elseif (args.rtype == 'it')
        fm = zeros(nim,nim,nslices);
        stsl = 0; ensl = nslices; 
        fm_penord = 2; fm_log2beta = -4.5; fm_niter = 300; fm_reps = 2;
        for slnum = stsl:ensl-1
            for ii = 1:fm_reps
                for phnum = 0:1
                    % do field map recon here
                    dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
                    immap = zeros([nim nim]);
                    for coilnum = 1:ncoils
                        % do standard recon
                        doiterrecio2;
                         
                        if (phnum == 0)
                            imkbsv = imkb;
                            if (scaninfo.revflg == 2)
                                if (im1nz == 1)
                                    im1sv = im1;
                                end
                            end
                        else
                            imkb = imkb .* conj(imkbsv);
                            if (scaninfo.revflg == 2)
                                if (im1nz == 1)
                                    imkb = imkb + im1 .* conj(im1sv);
                                end
                            end
                            immap = immap + imkb;
                        end % if (phnum == 0)
                    end % for coilnum
                end % for phnum
                mnfm = mean(immap(:));
                fmtmp = angle(immap.*exp(-i*angle(mnfm)));
                fmmag = sqrt(abs(immap));
                %fmsm = angle(conv2(fmmag.*exp(i*fmtmp),ones(args.fsize),'same'))+angle(mnfm);
                ticker reset;
                fmsm = mri_phase_denoise(fmmag.*exp(i*fmtmp), 'l2b', fm_log2beta, ...
                       'order', fm_penord, 'niter', fm_niter, 'wi_ml', true, 'pl', 0) + angle(mnfm);
                fm(:,:,slnum+1) = fmsm/(scaninfo.mapdel*1e-6)/2/pi; % in Hz
  
                % Clear the objects before redoing the field map estimation
                clear ito*;
            end % for field map estimate repetitions
            disp(sprintf('Field map for slice #%d done',slnum+1));
        end % for slnum
    end
    save fmfile fm;
else % do recon
    rtype = args.rtype;
    if args.usemap == 1
        load fmfile;
    else
        fm = zeros([nim nim nslices ncoils]);
    end
    if (iph == 0)
        stph = ((1-zflg)*mapdda*unfact);
        enph = nphases*unfact;
    else
        stph = iph-1+mapdda;
        enph = iph+mapdda;
    end
    if (isl == 0)
        stsl = 0; ensl = nslices;
    else
        stsl = isl-1; ensl = isl;
    end

    imar = zeros([nim nim nslices]);
    ahdr = make_ahdr(args,scaninfo);

    if args.allcoils == 1
        imarall = zeros([nim nim nslices ncoils]);
    end % if args.allcoils == 1
    if args.svraw == 1
        imarraw = zeros([OVER*nim OVER*nim nslices ncoils]);
    end % if args.svraw == 1
    for phnum = stph:enph-1
        phnumout = phnum+1-mapdda*unfact;
        for slnum = stsl:ensl-1
            % load data for for/rev, spirals, all coils
            dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
            for coilnum = 1:ncoils 
                % run recon
                if (args.rtype == 'cp') | (args.rtype == 'ft') % do ft or cp recon
                    docprec;
                end
                if (args.rtype == 'it') % do iter recon
                    doiterrecio2;
                end

                % combine if necessary
                if (scaninfo.revflg == 2)
                    if (im1nz == 1)
                        imkb = abs(abs(imkb) + i*abs(im1));
                    end
                end

                if args.flipx == 1
                    imkb = flipud(imkb);
                end
                if args.flipy == 1
                    imkb = fliplr(imkb);
                end
                if args.svraw == 1
                    imarraw(:,:,slnum+1,coilnum) = datgr;
                end % if args.svraw == 1
                if args.allcoils == 1
                    imarall(:,:,slnum+1,coilnum) = imkb;
                end % if args.allcoils == 1
                if ncoils > 1
                    % do L2 norm for multicoil
                    if coilnum == 1
                        imout = abs(imkb);
                    else
                        imout = abs(imout + i.*abs(imkb));
                    end
                    imkb = imout;
                end % if ncoils > 1
            end % for coilnum = 1:ncoils
            % now output indivual images image
            imar(:,:,slnum+1) = imkb;
            if strcmp(args.outfmt, 'sl')
                basestr = sprintf('sl%d.%03d', slnum, phnum);
                writeim(basestr, abs(imkb));
            end % if args.outfmt = 'sl'
        end % for slnum = stsl:ensl-1

        % write out full arrays
        basestr = sprintf('vol_e%d_%s_%04d',scaninfo.exam,datestr,phnum+args.imoffset);
        if args.svraw == 1
            eval(['save ' basestr '_raw imarraw;']);
        end % if args.svraw == 1
        if args.allcoils == 1
            eval(['save ' basestr '_all imarall;']);
        end % if args.allcoils == 1
        if strcmp(args.outfmt,'nifti')
            niidata(nii_start:nii_stop) = imar(:);
            nii_start = nii_stop + 1;
            nii_stop = nii_stop + pixperframe;
        end
        if strcmp(args.outfmt, 'anal')
            write_hdr([basestr '.hdr'],ahdr);
            write_img([basestr '.img'],abs(imar),ahdr);
        elseif strcmp(args.outfmt, 'mat')
            if args.complex == 0
                imar = abs(imar);
            end
            eval(['save ' basestr ' imar;']);
        end % if args.outfmt == 'anal'

    end %for phnum = stph:enph-1
end % if (args.map == 1)

fclose(fid);

if strcmp(args.outfmt,'nifti')
    nhdr = avw2nii_hdr(ahdr);
    nhdr.dim(1) = 4;
    nhdr.dim(5) = nphases -1;
    nhdr.dim(5) = nphases - mapdda;  % LHG - if you do a field map, the first one doesn't get counted.
    write_nii(sprintf('%s.nii',basestr), abs(niidata), nhdr,0);

    if args.complex==1
        write_nii(sprintf('p_%s.nii',basestr), 1000*angle(niidata), nhdr,0);
    end
end

