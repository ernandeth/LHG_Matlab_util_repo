function kdephaser(pfilename ,  showData )
% function kdephaser(pfilename , showData )
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
%  nslices=info1(3);
%  nechoes=info1(4);
%  nexcitations=info1(5);
%  nframes=info1(6);
%  framesize=info1(9);
%  psize=info1(10);
%

RAWHEADERSIZE=61464;
HEADERSIZE=61464;
endian = 'ieee-le';
RDBRAWHEADEROFFSET=30;  % For version 4 use 30.
LOC1 = RDBRAWHEADEROFFSET + 34;
LOC2 = RDBRAWHEADEROFFSET + 70;
LOC3 = RDBRAWHEADEROFFSET + 186;

% reading the header info in the P file:
if isstr(pfilename),
	fid=fopen(pfilename,'r', endian);
	if (fid<3), disp(['Could not open file! ...']); end;
else
	fid=pfilename;
end;

hdr=fread(fid, HEADERSIZE, 'uint8');

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

nfidsperslice = nechoes * (nframes+1)-1;
nfids = nslices * nechoes * (nframes+1)-1;
BASELINESIZE = framesize

if rev_for==1
	framesize=framesize*2;
	nfidsperslice = nfidsperslice/2;
	nfids = nfids/2;
end

% data type here
if (psize==4),
	fmt=['long'];
else,
	fmt=['short'];
end;
if (nargin==1),
	info1
	info2
	info3
end

% open the output file and stick the header in it
ofp=fopen(sprintf('f_%s', pfilename) ,'wb', endian);
fwrite(ofp, hdr, 'uint8');

Nspikes=0;
spikesperslice = [];

% move file pointer to where data begins
status=fseek(fid,RAWHEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;

maxes = zeros(1,nslices);
% do one slice at a time
for sl=1:nslices

	% read the complex data in
	baseline = fread(fid, 2*BASELINESIZE, fmt);
	raw = fread(fid, [2*framesize nfidsperslice], fmt);
	out = raw;


	%if ~isempty(find(kslices==sl))

	fprintf('\rslice: %d of %d ', sl, nslices);
	% turn into complex data

	comraw = complex(raw(1:2:end, :), raw(2:2:end, :) );
	magraw = abs(comraw);
	phsraw = angle(comraw);
	out2 = comraw;

	maxes(sl) = max(magraw(1,4));

	% find the difference in phase at the center of k-space and use that as
	% a correction.
	Npts = 5;
	
	for n=1:2:nfidsperslice
		phasecorr1 = angle( mean(comraw(1:Npts,n)) / mean(comraw(1:Npts,1)) );
		phasecorr2 = angle( mean(comraw(1:Npts,n+1)) / mean(comraw(1:Npts,2)) );
		
		out2(:,n) = comraw(:,n) .*exp(i * phasecorr1);
		out2(:,n+1) = comraw(:,n+1) .*exp(i * phasecorr2);
	end

	if (showData==1),
		subplot(221), imagesc(magraw); title('input magnitude');
		subplot(222), imagesc(phsraw); title('input phase')
		subplot(224), imagesc(angle(out2)); title('output phase')
		subplot(223), title(['slice ' num2str(sl) ' input' ])
		drawnow
		%colorbar

	end
	
	% split back into real and imaginary
	out(1:2:end,:) = real(out2);
	out(2:2:end,:) = imag(out2);

	out = out';
	% stick the baseline in before the data for the slice
	fwrite(ofp, baseline, fmt);
	for ct=1:size(out,1)
		fwrite(ofp, out(ct,:), fmt);
	end

end

plot(maxes);

if isstr(pfilename),
	fclose(fid);
end;
fclose(ofp);


