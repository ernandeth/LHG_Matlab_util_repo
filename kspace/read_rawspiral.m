function raw_timeseries = read_rawspiral(pfilename ,  showData )
% function raw_timeseries = read_rawspiral(pfilename , showData )
%
% P file is organized as follows:
%  coil 1:
%  slice 1: baseline interleaves time1, interleaves time 2, etc...
%  slice 2: baseline interleaves time1, interleaves time 2, etc...
%  etc...
%  coil 2: same as above but some random skip between them
%  if #leaves and #points is odd then there is an extra baseline at end
%

[args,scaninfo,kinfo] = rec_setup1(pfilename);

HEADERSIZE=scaninfo.headersize;
endian = 'ieee-le';

% reading the header info in the P file:
if isstr(pfilename),
	fid=fopen(pfilename,'r', endian);
	if (fid<3), disp(['Could not open file! ...']); end;
else
	fid=pfilename;
end;

nslices = scaninfo.nslices;
nechoes = scaninfo.npr;
nframes = scaninfo.nphases;
framesize = scaninfo.ndat;

nfidsperslice = nechoes * (nframes+1)-1;
nfids = nslices * nechoes * (nframes+1)-1;

nfidsperslice = nechoes * nframes;
nfids = nslices * nechoes * nframes;

fmt=['short'];

% allocate space:
raw_timeseries = zeros(nfidsperslice, nslices*framesize);

% move file pointer to where data begins
status=fseek(fid,HEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;

for sl=1:nslices

	% read the complex data in
	baseline = fread(fid, 2*framesize, fmt);
    
    % read one slice (all the frames)
	raw = fread(fid, [2*framesize nfidsperslice], fmt);
	
	fprintf('\rslice: %d of %d ', sl, nslices);
	% turn into complex data

	comraw = complex(raw(1:2:end, :), raw(2:2:end, :) );

    raw_timeseries( :, framesize*(sl-1)+1: framesize*sl) = comraw';

    if showData
        imagesc(abs(comraw')); drawnow
    end
end

fclose(fid);



