function result = orthoq(bname, image, windfact)
%function result = orthoq(bname, windfact)

bxname= sprintf('%sx.img',bname);
byname= sprintf('%sy.img',bname);
bzname= sprintf('%sz.img',bname);

bx = read_img2(bxname);
by = read_img2(byname);
bz = read_img2(bzname);
c = read_img2(image);

UR = 1;
% undersample:
bx = bx(1:UR:end, 1:UR:end, 1:UR:end);
by = by(1:UR:end, 1:UR:end, 1:UR:end);
bz = bz(1:UR:end, 1:UR:end, 1:UR:end);


b = sqrt(bx.^2 + by.^2 + bz.^2);
a = b;
if exist('windfact') == 0, 
  amin = min(a(:));
  amax = max(a(:));
  minmax = [amin,amax];
  a = (a  - amin);
else
  amin = windfact(1);
  amax = windfact(2);
  minmax = [amin,amax];
  a = (a  - amin);
  %a = a .* (a > 0);
end
a = (a)./(amax-amin).*64;
mask=ones(size(b));
mask(find(b < amin )) = 0;
mask(find(b > amax )) = 0;
b = mask .*a;
bx = mask.*bx;
by = mask.*by;
bz = mask.*bz;

x = 2;y=2; z=2;
colormap gray
[fig1, fig2, fig3] = oq([],c, bx,by,bz,x,y,z,0);
subplot(224), hist(b(:),100) 
button=1;
while (button ~= 122)  % until they press Z.  122 is ASCII for the Z key


    [i j button] = ginput(1);
    i=round(i);j=round(j);
    fig = floor(gca);
    switch(fig)
        case floor(fig1)
            x=j;
            y=i;
        case floor(fig2)
            z=j;
            x=i;
        case floor(fig3)
            y=i;
            z=j;
    end
    [fig1, fig2, fig3] = oq([],c,bx,by,bz,x,y,z,0);

end
