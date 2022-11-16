function indices=despiker_interp(pfilename , spike_threshold, fastrecflag, year, skip_frame)
% Usage ... indices = despiker(pfile, spike_threshold, fastrecflag, [,year] [,skip_frame])
% 
%  modified by W. Yau 6/27/2008 
%   in
%       start_frame     # of frames to skip from front (defaults to 0, i.e. no frames skipped)
%
% search for spikes in the kspace data and replace them by the mean of
% that point in kspace over the time series.
%
% a spike is defined as a point that has a value spike_threshold*stdDev greater than the
% mean of that spike location.
%
% P file is organized as follows:
%  coil 1:
%  slice 1: baseline interleaves time1, interleaves time 2, etc...
%  slice 2: baseline interleaves time1, interleaves time 2, etc...
%  etc...
%  coil 2: same as above but some random skip between them
%  if #leaves and #points is odd then there is an extra baseline at end
%
%
%  Slice 1 Frame 0 is baseline (example, expect more)
%  Note: there is frame 0 but not slice 0!
%
%  For sf3d - slc is the phase encode of interest and frm
%  would be the spiral shot of interest (0 is the baseline).
%
%  nslices=info1(3);
%  nechoes=info1(4);
%  nexcitations=info1(5);
%  nframes=info1(6);
%  framesize=info1(9);
%  psize=info1(10);
%  

% Check if scanner software version specified
fid = fopen(pfilename,'r');
s_hdr = read_gehdr(fid);
rev = s_hdr.rdb.rdbm_rev; 

switch round(rev)
    case 9
        % default if year is not specified....
        RAWHEADERSIZE=61464;
        endian = 'ieee-le';
        if nargin>4
            switch year
                case 2003
                    RAWHEADERSIZE=39984;
                    endian = 'ieee-be';
                case 2004
                    RAWHEADERSIZE=60464;
                    endian = 'ieee-be';
                case 2005
                    RAWHEADERSIZE=61464;
                    
            end
            if (year>=2005)
                RAWHEADERSIZE=61464;
            end
        end
        
    case 20
        RAWHEADERSIZE = 149788;
        endian = 'ieee-le';
    otherwise
        error('fmrilab:despiker:invalidrev','Invalid format rev %d',rev);
end
HEADERSIZE = RAWHEADERSIZE;

if ischar(pfilename),
    fid=fopen(pfilename,'r', endian);
    if (fid<3), disp(['Could not open file! ...']); end;
else
    fid=pfilename;
end;

RDBRAWHEADEROFFSET=30;  % For version 4 use 30.
LOC1 = RDBRAWHEADEROFFSET + 34;
LOC2 = RDBRAWHEADEROFFSET + 70;
LOC3 = RDBRAWHEADEROFFSET + 186;

% Get header information
hdr=fread(fid, HEADERSIZE, 'uint8');
switch round(rev)
    case 9
        status = fseek(fid,LOC1,'bof');
        if (status), disp(['Could not seek to file location 1! ...']); end;
        info1 = fread(fid,10,'short');
        
        status = fseek(fid,LOC2,'bof');
        if (status), disp(['Could not seek to file location 2! ...']); end;
        info2 = fread(fid,6,'short');
        
        status = fseek(fid,LOC3,'bof');
        if (status), disp(['Could not seek to file location 3! ...']); end;
        info3 = fread(fid,50,'float');
        nShots = info3(5);
        rev_for =info3(12)
        
        nslices = info1(3);
        nechoes = info1(4);
        nexcitations = info1(5);
        nframes = info1(6);
        framesize = info1(9);
        psize = info1(10);
        
        % Get number of coils (should be 1)
        [args,info] = rec_setup1(pfilename,'m','l','fy','n',64);
        ncoils = info.ncoils;
    case 20
        %rev_for = ??  [in v9 is info.s1.user11]
        rev_for = s_hdr.rdb.user11;
        nslices = s_hdr.rdb.nslices;
        nechoes = s_hdr.rdb.nechoes; % s_hdr.image.numecho;
        nframes = s_hdr.rdb.nframes;
        framesize = s_hdr.rdb.frame_size;
        psize = s_hdr.rdb.point_size;
        
        % Account for number of coils
        info = useoriginalfieldnames(s_hdr);
        ncoils = info.s1.dab(2)-info.s1.dab(1)+1;
end
nfidsperslice = nechoes * (nframes+1)-1;
nfids = nslices * nechoes * (nframes+1)-1;
baselinesize = framesize

if rev_for==1
    framesize=framesize*2;
    nfidsperslice = nfidsperslice/2;
    nfids = nfids/2;
end

status=fseek(fid,RAWHEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;

% data type here
if (psize==4),
    fmt=['long'];
else,
    fmt=['short'];
end;

indices=[0:nfidsperslice-1; zeros(1,nfidsperslice)];

% open the output file and stick the header and baseline in it
ofp=fopen(sprintf('f_%s', pfilename) ,'wb', endian);
fwrite(ofp, hdr, 'uint8');

% Do one coil at a time
for c = 1:ncoils
    
    % do one slice at a time
    for sl=1:nslices
        
        % read the complex data in
        baseline=fread(fid, 2*baselinesize, fmt);
        
        raw = fread(fid, [2*framesize nfidsperslice], fmt);
        raw = raw';
        re_raw = raw(:, 1:2:2*framesize);
        im_raw = raw(:, 2:2:2*framesize);
        raw = complex(re_raw, im_raw);
        
        
        if (fastrecflag),
            %demodulation?
            th = 2*pi*[1:framesize]'/4;
            th = kron(th , ones(1,size(raw,1)));
            draw = raw.*exp(-i*th');
        else
            draw =raw;
        end;
        
        out = zeros(size(draw));
        
        if skip_frame > 0
            out = draw((1+skip_frame):end,:);
        else
            out = draw;
        end
        
        hbadrows = zeros(1, size(out,1));
        % now do the filter (very crude stuff)
        
        % filter 1: amplitudes (valur version)
        aout = abs(out);
        meansa = mean(aout,1);
        stds2 = std(out,0,1);
        
        for col=1:size(out,2)
            badrows = find(aout(:,col) > ...
                (meansa(col) + spike_threshold*stds2(col)) );
            
            goodrows = find(aout(:,col) <= ...
                (meansa(col) + spike_threshold*stds2(col)) );
            
            
            lastrow = length(badrows);
            if (badrows) > 0
                if badrows(1) == 1
                    out(1,col) =  interp1(goodrows, out(goodrows,col), 1, 'nearest', 'extrap');
                end
                
                if badrows(lastrow) == size(out,1)
                    out(size(out,1),col) =  interp1(goodrows, out(goodrows,col), size(out,1), 'nearest', 'extrap');
                end
                
                out(:,col) = interp1(goodrows, out(goodrows,col), 1:size(out,1), 'linear', 'extrap');
                hbadrows = hbadrows + hist(badrows, 1:size(out,1));
                
                if badrows(1) == 1
                    out(1,col) =  interp1(goodrows, out(goodrows,col), 1, 'nearest', 'extrap');
                end
                
                if badrows(lastrow) == size(out,1)
                    out(size(out,1),col) =  interp1(goodrows, out(goodrows,col), size(out,1), 'nearest', 'extrap');
                end
            end
            
        end
        
        if skip_frame > 0
            indices(2,:) = indices(2,:) + [zeros(1,skip_frame) hbadrows];
            out = [draw(1:skip_frame,:); out];  % add back the zeroth frame of this slice
        else
            indices(2,:) = indices(2,:) + hbadrows;
        end
        
        if (fastrecflag),
            % modulation?
            out = out.*exp(i*th');
        end;
        
        re_out = real(out);
        im_out = imag(out);
        
        % put the data in the right format for writing
        out2=zeros(nfidsperslice , framesize*2);
        for count=1:framesize
            out2(:,2*count-1) = re_out(:,count);
            out2(:,2*count) = im_out(:,count);
        end
        
        fwrite(ofp, baseline, fmt);
        for ct=1:size(out2,1)
            %fprintf('\rWriting point %d of %d ..slice %d',ct,size(out2,2),sl);
            
            fwrite(ofp, out2(ct,:), fmt);
        end
        
    end
    
end

if isstr(pfilename),
    fclose(fid);
end;
fclose(ofp);


