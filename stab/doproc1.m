stnum = 63;  % study number
slnum = 7;   % slice number
stser = 3;   % series number
xloc = 32;   % roi loc 
yloc = 29;   % roi loc 

nsl = 12;
offset = nsl + slnum + 1;
header = 7904;
numim = 600;
xsz = 3;
ysz = 3;


im = 1;
imnum = im*nsl + slnum;
sernum = floor(imnum/999)*20 + stser;
subim = imnum - floor(imnum/999)*999 + 1;
str = sprintf('%05d/%03d/I.%03d',stnum,stser,subim);
b = readimsw(str,[64 64],header);
a1 = ones(size(b));
a1(xloc:xloc+xsz-1,yloc:yloc+ysz-1) = 1.2.*ones([xsz ysz]);
show((b.*a1)');
str = sprintf('study %05d, series %03d, slice %d (%dx%d)',stnum,sernum,slnum,xsz,ysz);
title(str);
pause
for im = 1:numim;
  imnum = im*nsl + slnum;
  sernum = floor(imnum/999)*20 + stser;
  subim = imnum - floor(imnum/999)*999 + 1;
  str = sprintf('%05d/%03d/I.%03d',stnum,sernum,subim);
  a = readimsw(str,[64 64],header);
  tc(im) = sum(sum(a(xloc:xloc+xsz,yloc:yloc+ysz)));
  ser(im) = sernum;
end
pe = std(tc)/mean(tc)*100
pe2 = std(detrend(tc))/mean(tc)*100
plot(tc./mean(tc))
str = sprintf('study %05d, series %03d, slice %d, err (%dx%d) = %.2f',stnum,stser,slnum,xsz,ysz,pe);
title(str);
ylabel('normalized image intensity')
xlabel('image number/time (s)');
