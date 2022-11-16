function y=readim(a,f,h)
% usage .. readim(a,f,h);
% displays matrix "a" is a string that names the input file
% and "f" is the optional size - default is [128 128]
% and "g" is the optional header size - default is 0 bytes

if exist('f') == 0, 
  f = [64  64 ];
end
if exist('h') == 0, 
  h = 0;
end

fid = fopen(a,'r');
if fid == -1
y = [];
else
tmp = fread(fid,h,'int8');
y = fread(fid,f,'short');
fclose(fid);
end
