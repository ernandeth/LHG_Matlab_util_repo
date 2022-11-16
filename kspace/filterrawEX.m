function [out,re_out,im_out,info1,info2,info3,ind]=filterraw(pfilename, bad_range ,fastrecflag, showData)
% Usage ... [out,qm,im,info1,info2,info3]=filterraw(pfile, bad_range,fastrecflag, showData)
%
% this file zeroes the frequencies indicated in the range (a vector os
% numbers)
% in all the frames of the indicated Pfile.
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

if(nargin<2)
	tmp= abs(fft(fftshift(readraw(pfilename,1)-readraw(pfilename,2))))+...
		abs(fft(fftshift(readraw(pfilename,3)-readraw(pfilename,4))))+...
		abs(fft(fftshift(readraw(pfilename,5)-readraw(pfilename,6))))+...
		abs(fft(fftshift(readraw(pfilename,7)-readraw(pfilename,8))));
	figure,plot(tmp);
	maxpoint=find(tmp==max(tmp(:)))
	bad_range=[maxpoint-2:maxpoint+2]
	drawnow, pause
end


RAWHEADERSIZE=60464;
RAWHEADERSIZE=61464;
% old P files :% RAWHEADERSIZE=40080;


HEADERSIZE=60464;
HEADERSIZE=61464;
% odl Pfiles:  HEADERSIZE=40080;

RDBRAWHEADEROFFSET=30;  % For version 4 use 30.
LOC1 = RDBRAWHEADEROFFSET + 34;
LOC2 = RDBRAWHEADEROFFSET + 70;
LOC3 = RDBRAWHEADEROFFSET + 186;
%OPXRESLOC = 74;
%OPYRESLOC = 76;
%NFRAMESLOC = 46;
%NSLICES =
%FRAMESIZELOC = 52;
%POINTSIZELOC = 54;

if isstr(pfilename),
	fid=fopen(pfilename,'r','l');
	if (fid<3), disp(['Could not open file! ...']); end;
else
	fid=pfilename;
end;

hdr=fread(fid, HEADERSIZE, 'char');

status=fseek(fid,LOC1,'bof');
if (status), disp(['Could not seek to file location 1! ...']); end;
info1=fread(fid,10,'short');

status=fseek(fid,LOC2,'bof');
if (status), disp(['Could not seek to file location 2! ...']); end;
info2=fread(fid,6,'short');

status=fseek(fid,LOC3,'bof');
if (status), disp(['Could not seek to file location 3! ...']); end;
info3=fread(fid,50,'float');

% if (nargin==1),
%
%     if (nargout==1), f=[info1;info2;info3]; else, f=info1; end;
%     qm=info2;
%     im=info3;
%     if (nargout==0), [info1;info2;info3], clear, end;
%     return
%
%else,

nslices=info1(3);
nechoes=info1(4);
nexcitations=info1(5);
nframes=info1(6);
framesize=info1(9);
psize=info1(10);

nfids = nslices * nechoes * (nframes+1)-1;

%    status=fseek(fid,-nfids*framesize*2*psize - 2*psize*48 ,'eof');
%    status=fseek(fid,-nfids*framesize*2*psize ,'eof');
%    status=fseek(fid,RAWHEADERSIZE + 2*framesize*psize,'bof');
status=fseek(fid,RAWHEADERSIZE ,'bof');
if (status),error(['Could not seek to file fid location! ...']); end;

% data type here
if (psize==4),
	fmt=['long'];
else,
	fmt=['short'];
end;


% read the complex data in
baseline=fread(fid, 2*framesize, fmt);
%baseline=fread(fid, 2*framesize, fmt);
%ind=ftell(fid)

raw = fread(fid, [2*framesize nfids/10], fmt);
raw = raw';
re_raw = raw(:, 1:2:2*framesize);
im_raw = raw(:, 2:2:2*framesize);
raw = complex(re_raw, im_raw);


%end;

if isstr(pfilename),
	fclose(fid);
end;

if nargin<3, fastrecflag=0; end;
if nargin<4, showData=0; end;

if (fastrecflag),
	%demodulation?
	th = 2*pi*[1:framesize]'/4;
	th = kron(th , ones(1,size(raw,1)));
	draw = raw.*exp(-i*th');
else
	draw =raw;
end;

% now do the filter (very crude stuff)
fdraw = fft(draw,[],2);
fout = fdraw;
fout(:,bad_range)= 0;

out = ifft(fout,[],2);




if (showData==1),
	%if (nargout==0),
	if (nargin==1),
		info1
		info2
		info3
	else,
		subplot(221)
		imagesc(abs(fdraw));
		title('FFT along k axis of input')
		colorbar
		subplot(222)
		imagesc(abs(fout));
		title('FFT along k axis of output')
		colorbar
		subplot(223)
		imagesc(abs(raw - out));
		title('abs(input - output)')
		colorbar
		for ct=2:size(raw,1)
			fprintf('FID number: %d ', ct)
			subplot(223)
			plot(abs(fdraw(ct,:)));
			title(sprintf('FFT of input FID # %d',ct));
			axis tight
			
			
			subplot(224)
			plot(abs(fdraw(ct-1,:) - fdraw(ct,:)));
			hold on
			plot(abs(fft(out(ct-1,:)) - fft(out(ct,:))),'r');
			axis tight
			hold off
			title('FFT (thisFID - lastFID)');
			%pause
			drawnow
			%set(gca, 'XLim',[1000 3000])
		end;
	end
end

if (fastrecflag),
	% modulation?
	out = out.*exp(i*th');
end;
%%%%%%%%%%
%out=raw;
%%%%%%%%%%

re_out = real(out);
im_out = imag(out);

% put the data in the right format for writing
out2=zeros(nfids , framesize*2);
for count=1:framesize
	out2(:,2*count-1) = re_out(:,count);
	out2(:,2*count) = im_out(:,count);
end

% Write the output P file
fp=fopen(sprintf('f_%s', pfilename) ,'wb', 'ieee-le');
fwrite(fp, hdr, 'char');
fwrite(fp, baseline, fmt);
for ct=1:size(out2,1)
	fwrite(fp, out2(ct,:), fmt);
end
fclose(fp);


