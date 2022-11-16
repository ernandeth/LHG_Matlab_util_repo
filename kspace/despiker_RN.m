function spikesperslice =despiker_RN(pfilename , spike_threshold, ncoils, showData ,year)
% Usage ... spikesperslice = despiker_R(pfile, spike_threshold,, ncoils,  showData [,year])
%
% search for spikes in the kspace data and replace them by the mean of
% its neighbors over time .
%
% a spike is defined as a point that has a value spike_threshold*stdDev greater than the
% mean of that spike location.
%
% the output is the spikes locations in columns: [ slice, time, kspace ]

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

% default if year is not specified....
RAWHEADERSIZE=61464;
HEADERSIZE=61464;
endian = 'ieee-le';
if nargin==5
	switch year
		case 2003
			RAWHEADERSIZE=40080;
			HEADERSIZE=40080;
			endian = 'ieee-be';

		case 2004
			RAWHEADERSIZE=60464;
			HEADERSIZE=60464;
			endian = 'ieee-be';
		case 2005
			RAWHEADERSIZE=61464;
			HEADERSIZE=61464;
			endian = 'ieee-le';

	end
	if (year>=2005)
		RAWHEADERSIZE=61464;
		HEADERSIZE=61464;
		endian = 'ieee-le';
	end
end

RDBRAWHEADEROFFSET=30;  % For version 4 use 30.
LOC1 = RDBRAWHEADEROFFSET + 34;
LOC2 = RDBRAWHEADEROFFSET + 70;
LOC3 = RDBRAWHEADEROFFSET + 186;

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
%keyboard
% if (nargin==1),
%
%     if (nargout==1), f=[info1;info2;info3]; else, f=info1; end;
%     qm=info2;
%     im=info3;
%     if (nargout==0), [info1;info2;info3], clear, end;
%     return
%
%else,
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



% open the output file and stick the header in it
ofp=fopen(sprintf('f_%s', pfilename) ,'wb', endian);
fwrite(ofp, hdr, 'uint8');

Nspikes=0;
spikesperslice = [];
sdiffdata = [];
%ncoils=8;

for cl=1:ncoils

	% do one slice at a time
for sl=1:nslices
	fprintf('\rcoil: %d slice: %d of %d (%d spikes so far)',cl,  sl, nslices, Nspikes);
	% read the complex data in
	baseline=fread(fid, 2*baselinesize, fmt);

	raw = fread(fid, [2*framesize nfidsperslice], fmt);
	raw = raw';

	% allocate space for output:
	out = raw;

	% Note :  rows = time,   cols = kspace
	% we also avoid the first time point, as it usually contains the field
	for col=1:size(out,2)
		% identify spikes at each k-space location:
		kpoint = detrend(raw(:,col));
		means = mean(kpoint);
		stds2 = std(kpoint);
		thresh = abs(means) + abs(spike_threshold* stds2);
		badrows = find(abs(kpoint) > thresh);
		% plot(abs(kpoint))
		% line([1 1000], [thresh thresh]);
		% drawnow
		for r=1:length(badrows)
			% replace spikes by a mean of its neighbors (in time) -
			% (linear interpolation)
			% make sure we're not overwriting the field map (the first row)!
			if (badrows(r)>2 & badrows(r) < (nfidsperslice -1))
				row = badrows(r);

				% out(row,col) = means(col);

				out(row, col) = 0.5*(raw(row-1, col) + raw(row+1, col)) ;
				Nspikes = Nspikes+1;
				spikesperslice = [spikesperslice ; sl badrows(r) col];
			end
		end
	end


	if (showData==1),
		%   if (sl==9),
		%if (nargout==0),
		if (nargin==1),
			info1
			info2
			info3
		else,
			subplot(221)
			imagesc(abs(raw));
			title(['slice ' num2str(sl) ' input' ])
			colorbar
			subplot(222)
			imagesc(abs(out));
			title('output')
			colorbar
			subplot(223)
			imagesc(abs(raw - out),[0 1]);
			title('abs(input - output)')
			colorbar
			subplot(224)
			diffraw = abs(diff(raw,1));
			imagesc(diffraw,[0 100]);
			title('DIfferences in neighbors')
			colorbar
			colormap(gray)
			
			drawnow
			str = sprintf('spikes_sl_%02d', sl);
			print(gcf, '-depsc', str);
			%{
            for ct=2:size(raw,1)
                subplot(224)
                %plot(out2(ct,1:2:end));
                %title(sprintf('output FID number: %d',ct));
                %hold on
                %plot(out2(ct,2:2:end),'r');
                %hold off
                plot(abs(raw(ct-1,:) - raw(ct,:)));
                hold on
                plot(abs(out(ct-1,:) - out(ct,:)),'r');
                hold off
                title('(thisFID - lastFID)');
                legend('in', 'out')
                %axis([0 4000 0 100])

                %pause
                subplot(223)
                plot(abs(raw(ct,:)));P64512.7
                hold on
                plot(abs(out(ct,:)),'r');
                hold off
                title(sprintf('slice: %d , FID # %d',sl, ct));
                legend('in', 'out')
                %axis([0 4000 0 100])
                drawnow
                pause
            end;
			%}
		end
	end
	if (showData==2)
		diffdata = abs(diff(raw,1));
		diffdata = sum(diffdata(2:end-1,:), 2);
		sdiffdata = [sdiffdata diffdata-mean(diffdata)];
		plot(sdiffdata);
		drawnow
		print -dps sumofDiffs_overTime
	end





	% stick the baseline in before the data for the slice
	fwrite(ofp, baseline, fmt);
	for ct=1:size(out,1)
		%fprintf('\rWriting point %d of %d ..slice %d',ct,size(out2,2),sl);
		fwrite(ofp, out(ct,:), fmt);
	end
end

end

if isstr(pfilename),
	fclose(fid);
end;
fclose(ofp);
fprintf('\nFound %d spikes .  Done', Nspikes);


