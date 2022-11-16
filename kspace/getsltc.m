function f=getsltc(path,slno,tlocs,imlocs,imsize,detord,ext)
% Usage ... f=getsltc(path,slno,tlocs,imlocs,imsize,detord,ext)

if nargin<5,
  imsize=[64 64];
end;
if nargin<7,
  ext='';
end;
if nargin<6,
  detord=0;
end;

if length(tlocs)==1,
  tlocs(2)=tlocs(1);
  tlocs(1)=1;
end;

for m=tlocs(1):tlocs(2),
  a=getslim(path,slno,m,imsize,ext);
  for n=1:size(imlocs,1),
    f(m,n)=a(imlocs(n,1),imlocs(n,2));
  end;
end;

if (detord>0),
  f=tcdetrend(f,detord); 
end;

