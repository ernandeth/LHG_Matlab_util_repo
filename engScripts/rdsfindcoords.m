function [inds] = rdsfindcoords(xc, yc, xres, yres, Nslices);
%
%function [inds] = rdsfindcoords(xc, yc, xres, yres, Nslices);
%

Nrows = floor(sqrt(Nslices));
Nrows = 3;
Ncols = ceil(Nslices/Nrows);

myrow = ceil(yc/yres) ;
mycol = ceil(xc/xres) ;

x_roi = round(rem(abs(xc),xres));
y_roi = round(rem(abs(yc),yres));
z_roi = Ncols * (myrow-1) + mycol;

%[myrow mycol x_roi y_roi z_roi]


InitialROIset=1;

if x_roi < 2, x_roi=2; end
if y_roi < 2, y_roi=2; end
if z_roi < 1, z_roi=1; end
if x_roi > xres-2,  x_roi=xres-1; end
if y_roi > xres-2,  y_roi=yres-1; end
if z_roi > Nslices, z_roi=Nslices; end

%
[xx,yy,zz] = ndgrid([x_roi-1:x_roi+1], [y_roi-1 : y_roi+1], z_roi*[ 1 1 1]);
nlist = length(xx(:));
fx = reshape(xx,nlist,1);
fy = reshape(yy,nlist,1);
fz = reshape(zz,nlist,1);

inds = sub2ind( [xres, yres, Nslices], fx, fy, fz);

return
