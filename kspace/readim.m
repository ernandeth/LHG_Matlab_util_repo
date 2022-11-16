function y=readim(a,f,h,type,compflag)
% usage .. readim(a,f,h,type,comptype);
% displays matrix "a" is a string that names the input file
% and "f" is the optional size - default is [128 128]
% and "g" is the optional header size - default is 0 bytes

%if exist('f') == 0, 
%  f = [128 128];
%end
%if exist('h') == 0, 
%  h = 0;
%end
%if exist('type') == 0,
%  type = 'short';
%end;

comp=computer;
if nargin<5, compflag=2; end;

if nargin<4, type='short'; end;
if nargin<3, h=0; end;
if nargin<2, f=[64 64]; end;

if (compflag==1),
  fid = fopen(a,'r','b');
elseif (compflag==2),
  fid = fopen(a,'r','l');
else,
  fid = fopen(a,'r');
end;
if (fid<3),
  disp(sprintf('Error: Can not open file %s',a));
end;
%tmp = fread(fid,h,'int8');
y = fread(fid,f,type);
fclose(fid);
