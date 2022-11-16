function f=getslvol(path,volno,imsize)
% Usage ... f=getslvol(path,volno,imno,[imsize volsize])

if length(imsize)==3,
  i1=1;
  iend=imsize(3);
else,
  i1=imsize(3);
  iend=imsize(4);
end;

for m=i1:iend,
  if volno<1000,
    fname=sprintf('%ssl%d.%03d',path,m,volno);
  else,
    fname=sprintf('%ssl%d.%04d',path,m,volno);
  end;
  %disp(fname);
  f(:,:,m)=readim(fname,imsize(1:2));
end;

