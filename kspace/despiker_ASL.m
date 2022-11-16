function spikesperslice =despiker_ASL(pfilename , spike_threshold, showData ,year)
% Usage ... spikesperslice = despiker_ASL(pfile, spike_threshold,  showData [,year])
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
if nargin==4
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

% do one slice at a time
for sl=1:nslices
	fprintf('\rslice: %d of %d (%d spikes so far)', sl, nslices, Nspikes);
	% read the complex data in
	baseline=fread(fid, 2*baselinesize, fmt);

	raw = fread(fid, [2*framesize nfidsperslice], fmt);
	%imagesc(abs(raw)); drawnow;
    raw = raw';
    


	% allocate space for output:
	out = zeros(size(raw));

	% Note :  rows = time,   cols = kspace
	% we also avoid the first time point, as it usually contains the field
	% for ASL, we need to spilt the data:
	r1 = raw(1:2:end,:);
	r2 = raw(2:2:end,:);

	out1 = r1;
	out2 = r2;

	% first do the odds
	for col=1:size(out1,2)
		% identify spikes at each k-space location:
		% first we detrend the data, then we identify the spikes
		kpoint = detrend(r1(:,col));

% 		dk= diff(kpoint);
% 		dmeans = mean(dk);
% 		dstds2 = std(dk);
% 		dthresh = abs(dmeans) + abs(spike_threshold* dstds2);
% 		dbadrows = find(abs(dk) > dthresh);


		means = mean(kpoint(2:end-1));
		stds2 = std(kpoint(2:end-1));
		thresh = abs(means) + abs(spike_threshold* stds2);
		badrows = find(abs(kpoint) > thresh);

		% 		if sl==3
		% 			plot(abs(kpoint))
		% 			line([1 length(kpoint)], [thresh thresh]);
		% 			drawnow
		% 		end

        % now replace the spikes by the mean of their neighbors
        for r=1:length(badrows)
            % replace spikes by a mean of its neighbors (in time) -
            % (linear interpolation)
            % make sure we're not overwriting the field map (the first row)!
            if (badrows(r)>2 & badrows(r) < (nfidsperslice/2 -2))
                % if the value is above threshold, check to see that the
                % derivative is too
                % 				if ~isempty(dbadrows==badrows(r)) | ...
                % 						~isempty(dbadrows== (badrows(r)-1))


                row = badrows(r);

                % out(row,col) = means(col);

                out1(row, col) = 0.5*(r1(row-1, col) + r1(row+1, col)) ;
                Nspikes = Nspikes+1;
                spikesperslice = [spikesperslice ; sl badrows(r) col];
                % 				end
            end
        end
    end


	% then do the evens
    for col=1:size(out2,2)
        kpoint = detrend(r2(:,col));

        % 		dk= diff(kpoint);
        % 		dmeans = mean(dk);
        % 		dstds2 = std(dk);
        % 		dthresh = abs(dmeans) + abs(0.5*spike_threshold* dstds2);
        % 		dbadrows = find(abs(dk) > dthresh);


        means = mean(kpoint(2:end-1));
        stds2 = std(kpoint(2:end-1),0,1);
        thresh = abs(means) + abs(spike_threshold* stds2);
        badrows = find(abs(kpoint) > thresh);

        % 		if sl==3
        % 			plot(abs(kpoint))
        % 			line([1 length(kpoint)], [thresh thresh]);
        % 			drawnow
        % 		end

        % now replace the spikes by the mean of their neighbors
        for r=1:length(badrows)
            % replace spikes by a mean of its neighbors (in time) -
            % (linear interpolation)
            % make sure we're not overwriting the field map (the first row)!
            if (badrows(r)>1 & badrows(r) < (nfidsperslice/2))
                % if the value is above threshold, check to see that the
                % derivative is too
                % 				if ~isempty(dbadrows==badrows(r)) | ...
                % 						~isempty(dbadrows== (badrows(r)-1))


                row = badrows(r);

                % out(row,col) = means(col);

                out2(row, col) = 0.5*(r2(row-1, col) + r2(row+1, col)) ;
                Nspikes = Nspikes+1;
                spikesperslice = [spikesperslice ; sl badrows(r) col];
                % 				end
            end
        end
    end


	% now merge the data
	out(1:2:end,:) = out1;
	out(2:2:end,:) = out2;


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
			imagesc(abs(raw - out));
			title('abs(input - output)')
			colorbar
			colormap(jet)
			drawnow
			str = sprintf('spikes_sl_%02d', sl);
			print(gcf, '-depsc', str);
			%
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
				plot(abs(raw(ct,:)));
				hold on
				plot(abs(out(ct,:)),'r');
				hold off
				title(sprintf('slice: %d , FID # %d',sl, ct));
				legend('in', 'out')
				%axis([0 4000 0 100])
				drawnow
				pause
			end;
			%
		end
	end




	% stick the baseline in before the data for the slice
	fwrite(ofp, baseline, fmt);
	for ct=1:size(out,1)
		%fprintf('\rWriting point %d of %d ..slice %d',ct,size(out2,2),sl);
		fwrite(ofp, out(ct,:), fmt);
	end

end

if isstr(pfilename),
	fclose(fid);
end;
fclose(ofp);
fprintf('\nFound %d spikes .  Done', Nspikes);


