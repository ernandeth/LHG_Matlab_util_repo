stnum = 53;
nsl = 12;
slnum = 3;
stser = 5;
offset = nsl + slnum + 1;
header = 7904;
imoff = 0;
numim = 300;
xloc = 34;
yloc = 34;
xsz = 3;
ysz = 3;
xsz2 = round(xsz/2);
xsz3 = xsz-xsz2-1;
ysz2 = round(ysz/2);
ysz3 = ysz-ysz2-1;

asz = 20;
asz2 = round(asz/2);
asz3 = asz-asz2-1;


im = 1;
imnum = im*nsl + slnum;
sernum = floor(imnum/999)*20 + stser;
subim = imnum - floor(imnum/999)*999 + 1;
str = sprintf('%05d/%03d/I.%03d',stnum,stser,subim);
a = readimsw(str,[64 64],header);
a1 = ones(size(a));
a1(xloc-xsz2:xloc+xsz3,yloc-ysz2:yloc+ysz3) = 1.2.*ones([xsz ysz]);
show((a.*a1)');
str = sprintf('study %05d, series %03d, slice %d (%dx%d)',stnum,sernum,slnum,xsz,ysz);
title(str);
pause
imev = zeros([64 64]);
imod = zeros([64 64]);
numev = 0;
numod = 0;
for im = 1:numim;
  imnum = (im+imoff)*nsl + slnum;
  sernum = floor(imnum/999)*20 + stser;
  subim = imnum - floor(imnum/999)*999 + 1;
  str = sprintf('%05d/%03d/I.%03d',stnum,sernum,subim);
  a = readimsw(str,[64 64],header);
  tc(im) = sum(sum(a(xloc-xsz2:xloc+xsz3,yloc-ysz2:yloc+ysz3)));
  ser(im) = sernum;
  if rem(im,2) == 1
    imev = imev + a;
    numev = numev+1;
  else
    imod = imod + a;
    numod = numod+1;
  end
end
avdif = (imev/numev - imod/numod);
avim = (imev/numev + imod/numod);
corav = avdif(xloc-asz2:xloc+asz3,yloc-asz2:yloc+asz3);
spatstd = std(corav(:))*sqrt((numev+numod)/(xsz*ysz))./mean(mean(avim(xloc-asz2:xloc+asz3,yloc-asz2:yloc+asz3)))*100

pe = std(tc)/mean(tc)*100
pedet = std(detrend(tc))/mean(tc)*100
plot(tc./mean(tc))
str = sprintf('study %05d, series %03d, slice %d, err (%dx%d) = %.2f',stnum,stser,slnum,xsz,ysz,pe);
title(str);
ylabel('normalized image intensity')
xlabel('image number/time (s)');
