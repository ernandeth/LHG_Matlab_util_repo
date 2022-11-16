function f=getslim(path,slno,imno,imsize,ext)
% Usage ... f=getslim(path,slno,imno,imsize,ext)
%
% For images written natively in Sun/Unix computers use [64 64 1]
% For images written natively in PC/Linux computers use [64 64 2]

if nargin<5, 
  ext=''; 
else, 
  if (~isempty(ext)), ext=[ext,'.']; end; 
end;

if imno<1000,
  fname=sprintf('%ssl%d.%s%03d',path,slno,ext,imno);
else,
  fname=sprintf('%ssl%d.%s%04d',path,slno,ext,imno);
end;
%disp(fname);

if nargin<4,
  imsize=[64 64];
end;
if length(imsize)==3,
  f=readim(fname,imsize(1:2),0,'short',imsize(3));
elseif length(imsize)==4,
  f=readim(fname,imsize(1:2),0,'float',imsize(3));
else,
  f=readim(fname,imsize);
end;

