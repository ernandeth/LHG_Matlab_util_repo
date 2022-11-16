function [raw, lro, lpe, thk, total_tr, te, ti, ref_pos, ss, sli_val] = epirecon( dirpath )

% function [raw, lro, lpe, thk, total_tr, te, ti, ref_pos, ss, sli_val] = epirecon( dirpath )
% Matlab script to do reconstruction of EPI data
%
% INPUT :
% dirpath : epi_fid directory path
% or pick the the epi_fid by GUI if no input argument
%
% OUTPUT :
% raw : copmplex data of the reconstructed image
% or write abs value in sdt file format if no output arguement
%
%   epi_fid
%      |Spikes Removal
%      |DC correction
%      |Even echo reverse
%      |tabc
%      |RdOut FT
%      |Ref scan correction
%      |Nav correction
%      |Shift PhEnc
%      |Zero-filling and Normalize to zero gain
%      |PhEnc FT
%      |(Partial Fourier Recon)
%      |(Zero Filling+Hanning Filter)
%      V
%   complex or SDT(abs)
%
%
% Cecil Chern-Chyi Yen chy25@pitt.edu
% Fedora 8 + Matlab 7.5

% 1/17/'08 *change the file name epi -> Epi
% 1/18/'08 +read shiftecho and negphases
% 1/27/'08 *FFT scaling
% 3/6/'08 *change from fft to ifft
% 3/10/'08 *switch RO and PE dimension
% 3/12/'08 *change disp to fprintf
% 3/31/'08 *rearrage remove ref scan, shift ppe and 2nd ifft
% 4/7/'08  *speed up exist by checking file name only
% 5/22/'08 *switch RO and PE dimension for the ppeshift matrix
% 7/23/'08 *preload the epi_fid into memory which is faster and clean up the code
% 7/27/'08 +output lro, lpe, thk, total_tr
% 7/28/'08 +option to output sdt
% 7/29/'08 +output te, ti and remove ref_pos in te, ti
%          *change output of scal from %g to %f
% 8/1/'08 +automatic zeor-filling when writing to sdt
% 8/3/'08 *take absolute value after zero-filling
% 8/4/'08 *reuse parameters like raw and index numbers
% 8/5/'08 *Change ifft to fft and move write sdt before fft along
%         *if input is absolute value,scal can be 65535
% 8/6/'08 *change DC correction to subsctract mean of np/16 and nnv/8 points and change nv to nnv
%         *fix ref corr: negtive sign
%         -don't remove ref scan like epiShaper
% 8/7/'08 *change navigator corr algorithm
% 8/8/'08 *switch RO and PE back at end
% 10/1/'08 +normalize to zero gain
% 10/22/'08 *move normalize to zero gain to the end
% 11/11/08 *change ' to .'
% 11/25/'08 *rewrite tabc
% 11/26/'08 *fix a bug in navigator echo corr (change nseg to p)
% 11/28/'08 *try to do nav and ref corr before tabc, but it is no better
% 12/2/'08 *swap table to fix a weird phase problem
% 12/3/'08 *fix silly bugs in rev and fftshift
%          *move DC corr before shuffling
% 12/4/'08 +display sfrq and some other parameters
% 12/5/'08 *do not swap RO and PE
% 12/7/'08 +display pss
% 12/8/'08 +sort slices
% 12/11/'08 +check acqstatus
%           +check console and set endian
% 12/15/'08 +uigetdir
%           +count recon time
% 12/18/'08 *make POCS working with shiftecho (negphase<8)
% 12/21/'08 +record diary to current time
% 1/13/'09 *change fgetl to fgets and strtok to sscanf
%          *change nargin to nargout in write SDT section
%          +output sli_val
% 1/14/'09 +output ss
% 1/26/'09 *read %f instead of %d when reading procpar for some parameters
% 1/28/'09 *fix shift ppe
% 1/29/'09 *don't fftshift on PE direction
% 2/2/'09 +hamming filter for zero-filling
% 2/5/'09 *fftshift on PE direction & fix shift ppe again
% 2/10/'09 *read complex fid inside the loop
% 2/12/'09 *change check pss algorithm
% 2/27/'09 *change to ifft again
%          *remove ref scan
% 3/2/'09 +threshold=std*3
% 3/24/'09 +check arraydim==1
%          +always write four dimension data
% 3/25/'09 +flip lr and ud for ifft
% 4/10/'09 *Do not flip PE direction in image domian
% 5/15/'09 *need to fftshift before and after fft
% 6/2/'09 +check seqcon *check total_tr
% 6/3/'09 +fftshift on RO before last fft *reverse shift ppe direction
% 6/24/'09 +spikes removal
% 6/26/'09 *use total_tr to caculate skip
% 7/24/'09 *detect ref_pos out of size, don't check spikes on last image
% 8/1/'09 *change j to 1i
% 8/5/'09 *for 2D data, add siz(3) and siz(4) in Half-Fourier
% 8/28/'09 +check arraydim >5 for caculate TSNR
% 9/29/'09 *change name of the log
% 5/19/'10 +process dualecho
% 5/23/'10 *check dualecho -remove check total_tr
%
% TODO: linear ramp sampling re-griding

if nargin < 1
    dirpath = uigetdir('~/data/', 'Pick the fid Directory');
end

%diary(['epirecon_' datestr(now, 'yymmddHHMMSS') '.log'])

tic

dirpath = [strtok(dirpath,'.') '.fid' ];

if( exist([dirpath '/fid.orig'], 'file') )
    pathFid = [dirpath '/fid.orig'];
elseif( exist([dirpath '/fid'], 'file') )
    pathFid = [dirpath '/fid'];
else
    error('fid and fid.orig do not exist')
end

if( exist([dirpath '/procpar.orig'], 'file') )
    pathProcpar = [dirpath '/procpar.orig'];
elseif( exist([dirpath '/procpar'], 'file') )
    pathProcpar = [dirpath '/procpar'];
else
    error('procpar and procpar.orig do not exist')
end

fprintf('RECONSTRUCTING %s:\n', pathFid)

%%%%%%%%%%%%%% parameter reading %%%%%%%%%%%%%%
% procpar HEADER
%   name,   subtype,    basictype,  max,    min,    step,   Ggroup, Dgroup, protection, active, intptr
%   %s,     %d,         %d,         %f,     %f,     %f,     %d,     %d,     %d,         %d

fprintf('reading procpar...')

fp = fopen( pathProcpar, 'rt' );

while 1
    line = fgets(fp);
    if (line == -1)
        break;
    else
        name = strtok(line);
        if (strcmp(name, 'np'))  % Number of data point (P)
            line=(sscanf(fgets(fp),'%d'))';
            np=line(2:end);
        elseif (strcmp(name, 'nv')) % Number of phase encoding steps (P)
            line=(sscanf(fgets(fp),'%d'))';
            nv=line(2:end);
        elseif (strcmp(name, 'ns')) % Number of slice to be acquired (P)
            line=(sscanf(fgets(fp),'%d'))';
            ns=line(2:end);
        elseif (strcmp(name, 'nseg')) % Number of segmentation (P)
            line=(sscanf(fgets(fp),'%d'))';
            nseg=line(2:end);
        elseif (strcmp(name, 'ne')) % Number of echoes to be acquired (P)
            line=(sscanf(fgets(fp),'%d'))';
            ne=line(2:end);
        elseif (strcmp(name, 'ref_pos')) % Referance position
            line=(sscanf(fgets(fp),'%d'))';
            ref_pos=line(2:end);
        elseif (strcmp(name, 'arraydim')) % Dimension of experiment (P)
            line=(sscanf(fgets(fp),'%d'))';
            arraydim=line(2:end);
        elseif (strcmp(name, 'negphases')) % Number of negative portion of k-space
            line=(sscanf(fgets(fp),'%d'))';
            negphases=line(2:end);
        elseif (strcmp(name, 'gain')) % Receiver gain (P) (dB)
            line=(sscanf(fgets(fp),'%d'))';
            gain=line(2:end);
        elseif (strcmp(name, 'acqstatus')) % Acquisition status (P)
            line=(sscanf(fgets(fp),'%d'))';
            acqstatus=line(2:end);
        elseif (strcmp(name, 'sli_val')) % Stimulation index (P)
            line=(sscanf(fgets(fp),'%d'))';
            sli_val=line(2:end);
        elseif (strcmp(name, 'ss')) % Number of steady state (P)
            line=(sscanf(fgets(fp),'%d'))';
            ss=line(2:end);
        elseif (strcmp(name, 'lro')) % FOV size for readout axis (P) (cm)
            line=(sscanf(fgets(fp),'%f'))';
            lro=line(2:end);
        elseif (strcmp(name, 'lpe')) % FOV size for phase-encode axis (P) (cm)
            line=(sscanf(fgets(fp),'%f'))';
            lpe=line(2:end);
        elseif (strcmp(name, 'thk')) % Slice thinkness (P) (mm)
            line=(sscanf(fgets(fp),'%f'))';
            thk=line(2:end);
        elseif (strcmp(name, 'total_tr')) % Total TR including multi-segements and multi-slices (ms)
            line=(sscanf(fgets(fp),'%f'))';
            total_tr=line(2:end);
        elseif (strcmp(name, 'tr')) % repeatition time (s)
            line=(sscanf(fgets(fp),'%f'))';
            tr=line(2:end);
        elseif (strcmp(name, 'ppe')) % Position of image center on phase encode axis (P) (cm) <- ppe is not working in epi, so need post-process
            line=(sscanf(fgets(fp),'%f'))';
            ppe=line(2:end);
        elseif (strcmp(name, 'te')) % echo time (s)
            line=(sscanf(fgets(fp),'%f'))';
            te=line(2:end);
        elseif (strcmp(name, 'ti')) % inversion time
            line=(sscanf(fgets(fp),'%f'))';
            ti=line(2:end);
        elseif (strcmp(name, 'sfrq')) % Transmitter frequency of observe nucleus (P) (Hz)
            line=(sscanf(fgets(fp),'%f'))';
            sfrq=line(2:end);
        elseif (strcmp(name, 'pss')) % Slice position (P) (cm)
            line=(sscanf(fgets(fp),'%f'))';
            pss=line(2:end);
        elseif (strcmp(name, 'pro')) % pro (P) (cm)
            line=(sscanf(fgets(fp),'%f'))';
            pro=line(2:end);
        elseif (strcmp(name, 'navecho')) % Navigator echo correction (Y or N)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            nav = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'dp')) % Double precision (Y or N)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            dp = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'petable')) % Table of phase encode (string)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            petable = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'shiftecho')) % shiftecho (Y or N)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            shiftecho = strtok(strjust(parm(1:end-1),'left'),'""');
	    %shiftecho = 'n';
        elseif (strcmp(name, 'presig')) % Preamplifier signal level selection (P) (h or l) (high-> gain-30dB)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            presig = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'console')) % System console type (string)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            console = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'seqcon')) % Acquisition loop control (string)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            seqcon = strtok(strjust(parm(1:end-1),'left'),'""');
        elseif (strcmp(name, 'nf')) % Number of FIDs (P)
            line=(sscanf(fgets(fp),'%d'))';
            nf=line(2:end);
        elseif (strcmp(name, 'dualecho')) % dualecho (Y or N)
            line = fgets(fp);
            [cnt, parm] = strtok(line);
            dualecho = strtok(strjust(parm(1:end-1),'left'),'""');
        end
    end
end

fclose( fp );
fprintf('done\n');

if ~strcmp(seqcon,'nccnn')
    error('Loop control is not nccnn! seqcon=%s',seqcon)
end

if nf ~=nv*ns
    error('nf does not macth nv*ns! nf=%g',nf)
end

if strcmp(nav,'y')
    nnv = nv-nseg;
else
    nnv = nv;
end

if ne~=1
    error('epirecon can not read multiple-echos FID')
end

if acqstatus(1)~=101
    error('Experiment does not completed! Done #%g, Error #%g',acqstatus(1),acqstatus(2))
end

if strcmp(console, 'inova')
    endian = 'ieee-be';  % default byteorder is Sun big-endian. If use Intel little-endian, change to ieee-le; datatype: WORD='int16'; REAL ='float32'; LWORD='int32'; COMPLEX='float32'
else
    error('Console type is not inova! console: %s', console)
end

if arraydim<2
    error('arraydim<2!')
end

total_tr=10;
skip=ceil(1800*4/total_tr);
if (ref_pos<arraydim-1 && ref_pos>0)||(ref_pos==0 && ss>skip)
    error('Some of the FID does not reach steady state because ref_pos=%g, ss=%g!',ref_pos,ss)
end

if ref_pos>arraydim-1
    error('no referance scan! ref_pos=%g',ref_pos)
end

fprintf('dimension=%gx%gx%gx%g(%gmmx%gmmx%gmm)\n', np/2, nnv, ns, arraydim ,lro*10, lpe*10, thk);

fprintf('pro/ppe/pss= %gmm /%gmm /%g', pro*10, ppe*10, pss(1)*10);
if length(pss)>1
    fprintf(',%g', pss(2)*10);
end
if length(pss)>2
    fprintf(',%g', pss(3)*10);
end
if length(pss)>3
    fprintf(',%g', pss(4)*10);
end
if length(pss)>4
    fprintf(',%g', pss(5)*10);
end
if length(pss)>5
    fprintf(',%g', pss(6)*10);
end
fprintf('mm\n');

fprintf('sfrq=%9.6fHz,TR=%gms,TE=%gms\n', sfrq, total_tr, te(1)*1000);

if( exist([ '/Users/hernan/matlab/recon/varianlib/tablib/' petable ], 'file') )
    pathpe = [ '/Users/hernan/matlab/recon/varianlib/tablib/' petable ];
elseif( exist( [ './' petable ], 'file') )
    pathpe = [ './' petable ];
else
    error('table: %s does not exist in /Users/hernan/matlab/recon/varianlib/tablib/ or ./', petable)
end

diary([strtok(dirpath,'.') '.log'])

%%%%%%%%%%%%%% epi_fid reading %%%%%%%%%%%%%%
% For more detail about header of fid file, check load_fid.m by Maj Hedehus at Varian
%
% MAIN HEADER
%   nblocks,    ntraces,    np,     ebytes, tbytes, bbytes, vers_id,    status, nbheaders
%   %d32,       %d32,       %d32,   %d32,   %d32,   %d32,   %d16,       %d16,   %d32
% status = s_data, s_spec, s_32, s_float, s_complex, s_hyper
%
% BLOCK HEADER
%   scale,  bstatus,    index,  mode,   ctcount,    lpval,  rpval,  lvl,   tlt
%   %d16,   %d16,       %d16,   %d16,   %d32,       %f32,   %f32,   %f32,  %f32

fprintf('reading epi_fid...')

fid = fopen(pathFid,'r', endian);
if strcmp(dp,'n')
    fid_raw = fread(fid,inf,'int16');
    main_hdr = 16;
    block_hdr = 14;
else
    fid_raw = fread(fid,inf,'int32');
    main_hdr = 8;
    block_hdr = 7;
end
fclose(fid);

if exist('dualecho','var') && dualecho == 'y'
    raw = 1i*ones(np/2,nv*2,ns,arraydim);
    for m=1:arraydim
        for n=1:nseg
            for p=1:ns
                for q=1:(nv*2/nseg)
                    offset = main_hdr+block_hdr+np*ne*(q-1)+(np*ne*nv/nseg)*(p-1)+(np*ne*(nv/nseg)*ns)*(n-1)+(block_hdr+np*ne*nv*ns)*(m-1);
                    raw(:,q+(nv/nseg)*(n-1),p,m) = fid_raw((offset+1):2:(offset+np)) + 1i*fid_raw((offset+2):2:(offset+np));
                end
            end
        end
    end
else
    raw = 1i*ones(np/2,nv,ns,arraydim);
    for m=1:arraydim
        for n=1:nseg
            for p=1:ns
                for q=1:(nv/nseg)
                    offset = main_hdr+block_hdr+np*ne*(q-1)+(np*ne*nv/nseg)*(p-1)+(np*ne*(nv/nseg)*ns)*(n-1)+(block_hdr+np*ne*nv*ns)*(m-1);
                    raw(:,q+(nv/nseg)*(n-1),p,m) = fid_raw((offset+1):2:(offset+np)) + 1i*fid_raw((offset+2):2:(offset+np));
                end
            end
        end
    end
end


clear fid_raw;

fprintf('done\n')

%%%%%%%%%%%%%% Split Dual Echo %%%%%%%%%%%%%%



%%%%%%%%%%%%%% Spikes Removal %%%%%%%%%%%%%%
% find the spikes defined as abs(change)>(0.1mean+5std)

if 0
if arraydim>5
    fprintf('detecting spikes in FID...')
    
    if ss<skip
        skip=skip-ss;
    else
        skip=0;
    end
    if ref_pos==arraydim-1
        raw_mean=sum(abs(raw(:,:,:,skip+1:arraydim-1)),4)/(arraydim-skip-1);
        temp=0;
        for tt=skip+1:arraydim-1
            temp=temp+(abs(raw(:,:,:,tt))-raw_mean).^2;
        end
        raw_std=sqrt(temp./(arraydim-skip-1));
        for tt=skip+1:arraydim-2
            for slc=1:ns
                [spike_x, spike_y]=find((abs(raw(:,:,slc,tt))>(1.1*raw_mean(:,:,slc)+5*raw_std(:,:,slc)))+(abs(raw(:,:,slc,tt))<(0.9*raw_mean(:,:,slc)-5*raw_std(:,:,slc))));
                if ~isempty(spike_x)
                    for m=1:length(spike_x)
                        fprintf(' @(%g,%g,%g,%g)',spike_x(m),spike_y(m),slc,tt)
                        raw(spike_x(m),spike_y(m),slc,tt)=(raw(spike_x(m),spike_y(m),slc,tt-1)+raw(spike_x(m),spike_y(m),slc,tt+1))/2;
                    end
                end
            end
        end
    elseif ref_pos==0
        raw_mean=sum(abs(raw(:,:,:,skip+2:arraydim)),4)/(arraydim-skip-1);
        temp=0;
        for tt=skip+2:arraydim
            temp=temp+(abs(raw(:,:,:,tt))-raw_mean).^2;
        end
        raw_std=sqrt(temp./(arraydim-skip-1));
        for tt=skip+2:arraydim-1
            for slc=1:ns
                [spike_x, spike_y]=find((abs(raw(:,:,slc,tt))>(1.1*raw_mean(:,:,slc)+5*raw_std(:,:,slc)))+(abs(raw(:,:,slc,tt))<(0.9*raw_mean(:,:,slc)-5*raw_std(:,:,slc))));
                if ~isempty(spike_x)
                    for m=1:length(spike_x)
                        fprintf(' @(%g,%g,%g,%g)',spike_x(m),spike_y(m),slc,tt)
                        raw(spike_x(m),spike_y(m),slc,tt)=(raw(spike_x(m),spike_y(m),slc,tt-1)+raw(spike_x(m),spike_y(m),slc,tt+1))/2;
                    end
                end
            end
        end
    end
    fprintf('done\n')
    fprintf('TSNR at center of K-Space=%g for the first slice\n',raw_mean(np/4+1,nnv/2+1,1)/raw_std(np/4+1,nnv/2+1,1))
end
end

%%%%%%%%%%%%%% DC Correction %%%%%%%%%%%%%%
num=4;

fprintf('removing DC offset(num=%g)...',num)

if ( ref_pos+1 == arraydim )
    pos = arraydim-1;
else
    pos = arraydim;
end

for m=1:nseg;
    for n=1:ns
        raw(:,nv/nseg*(m-1)+1:nv/nseg*m,n,:) = raw(:,nv/nseg*(m-1)+1:nv/nseg*m,n,:)-mean(mean(raw(np/2-num+1:np/2,nv/nseg*m-num+1:nv/nseg*m,n,pos),1),2);
    end
end

fprintf('done\n')

%%%%%%%%%%%%%% even echo reverse %%%%%%%%%%%%%%

fprintf('flipping rows...')

rev=[];
for m = 1:nseg;
    rev=[rev nv/nseg*(m-1)+2:2:nv/nseg*m];
end

raw(:,rev,:,:) = raw(np/2:-1:1,rev,:,:);

fprintf('done\n')

%%%%%%%%%%%%%% Read phase encode tables %%%%%%%%%%%%%%

fprintf('converting data in table order to linear order...')
pathpe
inPtr = fopen(pathpe, 'rt',endian);
inStr = fgets(inPtr);
table = fscanf(inPtr,'%d');
fclose(inPtr);

[inStr,tabc] = sort(table);
if strcmp(nav,'y') % swap in PE direction before fft
    tabc = [tabc(nnv/2+1:nnv); tabc(1:nnv/2); tabc(end-nseg/2+1:end); tabc(end-nseg+1:end-nseg/2)];
else
    tabc = [tabc(nv/2+1:end); tabc(1:nv/2)];
end

raw = raw(:,tabc',:,:);

fprintf('done\n')

%%%%%%%%%%%%%% Discrete Fast Fourier Transform Along RO %%%%%%%%%%%%%%

fprintf('ifft along readout direction...')
raw = ifft(ifftshift(raw,1),np/2,1)*np/2;
fprintf('done\n')

%%%%%%%%%%%%%% reference scan correction %%%%%%%%%%%%%%

if (ref_pos > arraydim-1) || (ref_pos < 0)
    fprintf('no referance scan! ref_pos=%g\n', ref_pos)
else
    fprintf('reference scan correction(ref_pos=%g)...',ref_pos)
    angRef = angle(raw(:,:,:,ref_pos+1));
    raw(:,:,:,ref_pos+1)=[];
    arraydim=arraydim-1;
    for n=1:arraydim
        raw(:,:,:,n) = raw(:,:,:,n) .* exp( -1i*angRef );
    end
    fprintf('done\n')
    clear angRef
end

%%%%%%%%%%%%%% navigator echo correction %%%%%%%%%%%%%%
% Ehman, Felmlee; Radiology 173:255 (1989)
% Hu, Kim; MRM 31:495 (1994)
% Kim et al.; MRM 35:895 (1996)

if strcmp(nav,'y')
    fprintf('navigator echo correction...')
    raw_nav = raw(:,1:nnv,:,:);
    for m=1:arraydim
        for n=1:ns
            for p=2:nseg
                raw_nav(:,nnv/nseg*(p-1)+1:nnv/nseg*p,n,m) = raw(:,nnv/nseg*(p-1)+1:nnv/nseg*p,n,m).*(exp(-1i*angle(raw(:,nnv+p,n,m)./raw(:,nnv+1,n,m)))*ones(1,nnv/nseg));
            end
        end
    end
    raw = raw_nav;
    %raw=raw(:,1:(end-nseg),:,:);
    clear raw_nav
    fprintf('done\n')
end



%%%%%%%%%%%%%% Shift ppe %%%%%%%%%%%%%%
ppe=0;
%ppe=-lpe/2;

if ppe ~= 0
    fprintf('shifting ppe=%gmm...',ppe*10)
    peShift = ones(np/2,1)*exp(-1i*2*pi*ppe*([0:nnv/2-1,-nnv/2:-1])/lpe);
    for m=1:arraydim
        for n=1:ns
            raw(:,:,n,m) = raw(:,:,n,m).*peShift;
        end
    end
    fprintf('done\n')
    clear peShift
end

%%%%%%%%%%%%%% sort the slices %%%%%%%%%%%%%%

[inStr,shuf] = sort(pss);

if sum(pss~=inStr)
    fprintf('sorting slices:pss=')
    for m=1:size(pss,2)
        fprintf(' %g',inStr(m)*10)
    end
    raw = raw(:,:,shuf,:);
    fprintf('mm...done\n')
end

%%%%%%%%%%%%%% normalize to zero gain %%%%%%%%%%%%%%

if strcmp(presig,'h')
    gain = gain-30;
end
fprintf('normalizing to zero gain(gain=%g)...',gain)
raw = raw*10^(-gain/20);
fprintf('done\n')

%%%%%%%%%%%%%% Write to SDT WORD File %%%%%%%%%%%%%%

if ~nargout
    if ndims(raw) == 2
        siz = [size(raw), 1, 1 ];
    elseif ndims(raw) == 3
        siz = [size(raw), 1 ];
    else
        siz = size(raw);
    end
    
    zfRO = 2^nextpow2(siz(1));
    zfPE = 2^nextpow2(siz(2));
    %     zfRO = 128;
    %     zfPE = 128;
    
    if (siz(1) ~= zfRO ) || ( siz(2) ~= zfPE )
        fprintf('zero-filling to %gx%g and ifft along phase encoding direction...', zfRO, zfPE)
        raw = fft(raw,[],1)/np*2;
        
        hamm2d=(25-21*cos(2*pi*(0:siz(1)-1)'/siz(1)))/46*(25-21*cos(2*pi*(0:siz(2)-1)/siz(2)))/46;
        hamm2d=fftshift(hamm2d*siz(1)*siz(2)/sum(sum(hamm2d)));
        for p=1:siz(3)
            for q=1:siz(4)
                raw(:,:,p,q)=raw(:,:,p,q).*hamm2d;
            end
        end
        
        raw_zf = zeros(zfRO,zfPE,siz(3),siz(4));
        
        raw_zf(1:siz(1)/2,1:siz(2)/2,:,:) = raw(1:siz(1)/2,1:siz(2)/2,:,:);
        raw_zf(end-siz(1)/2+1:end,1:siz(2)/2,:,:) = raw(end-siz(1)/2+1:end,1:siz(2)/2,:,:);
        raw_zf(1:siz(1)/2,end-siz(2)/2+1:end,:,:) = raw(1:siz(1)/2,end-siz(2)/2+1:end,:,:);
        raw_zf(end-siz(1)/2+1:end,end-siz(2)/2+1:end,:,:) = raw(end-siz(1)/2+1:end,end-siz(2)/2+1:end,:,:);
        
        raw = ifftshift(ifftshift(ifft2(raw_zf),1),2)*zfPE*zfRO;
        %raw = ifftshift(ifft2(raw_zf),1)*zfPE*zfRO;
        clear raw_zf;
        siz(1)=zfRO;
        siz(2)=zfPE;
        
        fprintf('done\n')
    else
        fprintf('ifft along phase encoding direction...');
        raw = ifftshift(ifftshift(ifft(raw,[],2),1),2)*nnv;
        %raw = ifftshift(ifft(raw,[],2),1)*nnv;
        fprintf('done\n')
    end
    
    raw = abs(raw);
    raw = raw(end:-1:1,end:-1:1,:,:);
    %scal = 65535/max(max(max(max(raw)))); % Caculate float32 convert to int16 ratio
    
    % Write sdt file
    fprintf('writing %s.sdt...',strtok(dirpath,'.'))
    fid = fopen([strtok(dirpath,'.'), '.sdt'],'w', endian);
    %fwrite(fid, raw*scal, 'int16'); % dataType
    fwrite(fid, raw, 'float32'); % dataType
    fclose(fid);
    
    % Write spr file
    fid = fopen([strtok(dirpath,'.'), '.spr'],'wt', endian);
    fprintf(fid,'numDim: 4\n');
    fprintf(fid,'dim: %s\n', num2str(siz));
    fprintf(fid,'interval: %f %f %f %f\n', lro/siz(1)*10, lpe/siz(2)*10, thk, tr);
    %fprintf(fid,'dataType: %s\n', 'WORD'); % dataType
    %fprintf(fid,'Real2WordScale: %f\n', scal);
    fprintf(fid,'dataType: %s\n', 'REAL'); % dataType
    fprintf(fid,'endian: %s\n', endian); % byteorder
    fclose(fid);
    fprintf('done\n')
    
    clear raw
    
else
    %%%%%%%%%%%%%% POCS reconstruction for shiftecho/half-fourier %%%%%%%%%%%%%%
    if (shiftecho == 'y') && (nav == 'n')
        fprintf('POCS reconstruction for negphases=%g...',negphases)
        
        % 1 Create a truncated and zero padded k-space data(AAzf) from the shiftecho data(AAraw)
        AAraw = fftshift(fftshift(fft(raw,[],1),1),2)/np*2;
        %AAraw = fftshift(fft(raw,[],1),1)/np*2;
        
        new_siz=(nnv-negphases)*2;
        AAzf = zeros([size(AAraw,1),new_siz,size(AAraw,3),size(AAraw,4)]);
        cmplx_new = AAzf;
        cmplx_old = AAzf;
        AAzf(:,nnv-negphases*2+1:nnv,:,:) = AAraw(:,1:negphases*2,:,:);
        
        % 2 Generate a low freq. phase map(lowph)
        AAlowph = angle(ifft2(ifftshift(ifftshift(AAzf,1),2),np/2,new_siz));
        
        % 3 Use the aymmetric plus zero-filling data as our inital magnitude data
        AAzf(:,nnv+1:end,:,:)=AAraw(:,negphases*2+1:end,:,:);
        
        % 4 Generate an initial complex data(cmplx)
        AAcmplx = abs(ifft2(fftshift(fftshift(AAzf,1),2),np/2,new_siz)).*exp(AAlowph*1i);
        
        for tt=1:arraydim
            for slc=1:ns
                threshold=sqrt(var(reshape(AAzf(:,:,slc,tt),size(AAzf,1)*size(AAzf,2),1)))*3;
                for itr=1:50
                    % 5 Create an intermediate k-space data
                    AAzf(:,:,slc,tt) = fft2(AAcmplx(:,:,slc,tt),np/2,new_siz);
                    
                    % 6 Subsitute part of the intermediate k-space data into shiftecho data(raw_new)
                    AAzf(:,:,slc,tt)=ifftshift(ifftshift(AAzf(:,:,slc,tt),1),2);
                    AAzf(:,nnv-negphases*2+1:end,slc,tt)=AAraw(:,:,slc,tt);
                    AAzf(:,:,slc,tt)=fftshift(fftshift(AAzf(:,:,slc,tt),1),2);
                    
                    % 7 fft to create new complex data
                    cmplx_new(:,:,slc,tt) = ifft2(AAzf(:,:,slc,tt),np/2,new_siz);
                    
                    if sum(sum(abs(cmplx_new(:,:,slc,tt)-cmplx_old(:,:,slc,tt))))<threshold
                        AAzf(:,:,slc,tt)=ifftshift(ifftshift(AAzf(:,:,slc,tt),1),2);
                        break
                    end
                    cmplx_old(:,:,slc,tt)=cmplx_new(:,:,slc,tt);
                    
                    % 10 Generate an intermediate complex data(cmplx) and iterate
                    AAcmplx(:,:,slc,tt) = abs(cmplx_new(:,:,slc,tt)).*exp(AAlowph(:,:,slc,tt)*1i);
                end
            end
        end
        
        % 9 Smoothen via uu-point Hanning filter and output the result to AAcmplx
        uu = 3;
        for pp = 1:uu
            AAzf(:,nnv-negphases*2+pp,:,:) = (AAraw(:,pp,:,:)*(1+cos(pi+pi*(pp-1)/uu)) + AAzf(:,pp,:,:)*(1+cos(pi*(pp-1)/uu)))/2;
        end
        
        siz=size(AAzf);
        if length(siz)<4
            siz(4)=1;
        end
        if length(siz)<3
            siz(3)=1;
        end
        hamm2d=(25-21*cos(2*pi*(0:siz(1)-1)'/siz(1)))/46*(25-21*cos(2*pi*(0:siz(2)-1)/siz(2)))/46;
        for p=1:siz(3)
            for q=1:siz(4)
                AAzf(:,:,p,q)=AAzf(:,:,p,q).*hamm2d;
            end
        end
        
        raw = fftshift(fftshift(ifft2(fftshift(fftshift(AAzf,1),2),np/2,new_siz),1),2)*new_siz*np/2;
        %raw = fftshift(ifft2(fftshift(fftshift(AAzf,1),2),np/2,new_siz),1)*new_siz*np/2;
        
        clear AAlowph AAraw AAcmplx cmplx_new cmplx_old
        fprintf(' done\n')
    else
        %%%%%%%%%%%%%% Discrete Fast Fourier Transform Along PE %%%%%%%%%%%%%%
        fprintf('ifft along phase encoding direction...');
        raw = ifftshift(ifftshift(ifft(raw,nnv,2),1),2)*nnv;
        %raw = ifftshift(ifft(raw,[],2),1)*nnv;
        fprintf('done\n')
    end
    raw=raw(end:-1:1,end:-1:1,:,:);
end

fprintf('epirecon finished in %5.2f sec\n',toc);

diary off
%figure,subplot(221),imagesc(real(raw(:,:,1,1))), colormap('gray'), grid on
%subplot(222),imagesc(imag(raw(:,:,1,1))), colormap('gray'), grid on
%subplot(223),imagesc(abs(raw(:,:,1,1))), colormap('gray'), grid on
%subplot(224),imagesc(angle(raw(:,:,1,1))), colormap('gray'), grid on
