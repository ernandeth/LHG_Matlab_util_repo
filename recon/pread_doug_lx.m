function [newdata,fiddata,bsln]=pread_doug_lx(filename)
% Returns the raw data, parsed to follwing format:
% 	newdat[2 ndat npr nphases nslices]
% where the first index is 1 for real, two for imaginary
% optional fiddata returns size of each dimension

RAWHEADERSIZE=39984;


% Reading header of original pfile
[info1,info2,info3]=fidread2(filename);

ndat=info1(9);
nslices=info1(3)/info3(5);
npr=info3(5);
nph=1;
%nphmult=info3(11);
nphases=info1(6);
nframes=(npr*nphases+1)*nslices;

fiddata(1) = ndat;
fiddata(2) = npr;
fiddata(3) = nphases;
fiddata(4) = nslices;

% Opening files
file_id=fopen(filename,'r');
if file_id<3, error('Could not open file...'); end;

disp(['Reading raw data']);
status=fseek(file_id,RAWHEADERSIZE,'bof');
if status, error('Could not seek past header location...'); end;

size([2*ndat nframes])

[raw_data,rcnt]=fread(file_id,[2*ndat nframes],'short');
if rcnt~=nframes*2*ndat, error('Could not read all raw data...'); end;
fclose(file_id);

% Removing baselines

disp(['Removing baselines']);


new=zeros([2*ndat npr*nphases*nslices]);
 
for x=1:(nph*nslices),

	startx=(x-1)*(npr*nphases +1) + 1;
	startn=(x-1)*(npr*nphases) + 1;
	bsln(:,x)=raw_data(:,startx);
	new(:,startn:(startn+(npr*nphases -1))) = raw_data(:,(startx+1):(startx+(npr*nphases)));
end;

% Reshaping

disp(['Reshaping']);

newdata = reshape(new,[2 ndat npr nphases nslices]);

