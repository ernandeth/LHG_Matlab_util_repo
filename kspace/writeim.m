function y=writeim(a,b,type)
% usage .. writeim(a,b,type);
% writes array b into file "a" 
% output is number of bytes written
if nargin<3, type='short'; end;
fid = fopen(a,'w');
y = fwrite(fid,b,type);
fclose(fid);
