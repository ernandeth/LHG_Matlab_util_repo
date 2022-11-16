function [raw_timeseries, scaninfo] = read_raw_3d(pfilename ,  showData )
% function raw_timeseries = read_raw_3d(pfilename , showData )
%
% Luis Hernandez-Garcia @ UM (10/23/13)
%
% P file is organized as follows:
%   coil 1:
%       slice 1: 
%             baseline 
%             Frame 
%                 interleave, 
%  coil 2: same as above but some random skip between them
%  if #leaves and #points is odd then there is an extra baseline at end
%
% The output of this program will look like a 2D matrix of complex numbers for each channel:
% only 1 channel implemented right now.
%   rows = frames
%   colums = data in all FIDS concatenated horizontally. 
%       (slice1 (shot1,shot2...) - slice2 (shot1, shot2 ...) - 
%
% if getTraj, then we read in grad.txt and rotmats.txt and form a
% trajectory for each shot
%

[args,scaninfo,kinfo] = rec_setup1(pfilename);  % get header info using Doug's program

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
fidsize = scaninfo.ndat;
ncoils = scaninfo.ncoils;

nfidsperslice = nechoes * nframes;
nfids = nslices * nechoes * nframes;

fprintf('\nN echoes: %d', nechoes);
fprintf('\nN frames: %d', nframes);
fprintf('\nN slices: %d', nslices);
fprintf('\nN coils: %d',  ncoils);

fprintf('\n\n');

fmt=['short'];

% allocate space:
raw_timeseries = zeros(nframes, nslices*nechoes*fidsize, ncoils);

% move file pointer to where data begins
status=fseek(fid,HEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;


for nc=1:ncoils
    for sl=1:nslices
        
        % read the complex data in
        baseline = fread(fid, 2*fidsize, fmt);
        
        % read one slice (all the frames)
        raw = fread(fid, [2*fidsize nfidsperslice], fmt);
        
        %fprintf('\rcoil: %d  slice: %d of %d ', nc, sl, nslices);
        % turn into complex data
        
        comraw = complex(raw(1:2:end, :), raw(2:2:end, :) );
        comraw2 = reshape(comraw, nechoes*fidsize, nframes)';
        
        raw_timeseries( :, fidsize*nechoes*(sl-1)+1: nechoes*fidsize*sl, nc) = comraw2;
        
        if showData
            imagesc(abs(comraw')); drawnow
            % plot(abs(comraw)); drawnow
        end
    end
end

fclose(fid);


return

    
