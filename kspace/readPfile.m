function raw = readPfile(pfilename, dims)
%

	
fid=fopen(pfilename,'rb','ieee-le');

fmt=['short'];
fmt = ['long'];
%HEADERSIZE = 61464;

fseek(fid, -4*dims(1)*dims(2)*dims(3)*dims(4),'eof');

raw = fread(fid, dims(1)*dims(2)*dims(3)*dims(4)*2, fmt);

re_raw = raw( 1:2:end);
im_raw = raw( 2:2:end);
raw = complex(re_raw, im_raw);

%raw = reshape(raw, dims(1), dims(2), dims(3));


if isstr(pfilename),
	fclose(fid);
end;





