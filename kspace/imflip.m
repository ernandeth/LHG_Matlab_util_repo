function y=imflip(im,flipdir)

if (nargin<2), flipdir=1; end;

if (ischar(flipdir)),
  if (flipdir=='x'), flipdir=1; 
  elseif (flipdir=='y'), flipdir=2; 
  else, flipdir=1;
  end;
end;

xdim=size(im,1);
ydim=size(im,2);

y=zeros(size(im));

if (flipdir==2),
  for m=1:xdim,
    for n=1:ydim,
      y(m,n)=im(m,ydim-n+1);
    end;
  end;
else,
  for n=1:ydim,
    for m=1:xdim,
      y(m,n)=im(xdim-m+1,n);
    end;
  end;
end;

