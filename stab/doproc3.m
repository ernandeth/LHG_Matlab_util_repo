stnum = 36;
nsl = 12;
stser = 2;
header = 7904;
numim = 600;
xloc = 31;
yloc = 29;
xsz = 3;
ysz = 3;

slnum = 1;
im = 1;
imnum = im*nsl + slnum-1;
sernum = floor(imnum/999)*20 + stser;
subim = imnum - floor(imnum/999)*999 + 1;
str = sprintf('%05d/%03d/I.%03d',stnum,stser,subim);
a = readim(str,[64 64],header);
a1 = ones(size(a));
a1(xloc:xloc+xsz-1,yloc:yloc+ysz-1) = 1.2.*ones([xsz ysz]);
show((a.*a1)');
str = sprintf('study %05d, series %03d, slice %d (%dx%d)',stnum,sernum,slnum,xsz,ysz);
title(str);
pause
for slnum=1:nsl
for im = 1:numim;
  imnum = im*nsl + slnum-1;
  sernum = floor(imnum/999)*20 + stser;
  subim = imnum - floor(imnum/999)*999 + 1;
  str = sprintf('%05d/%03d/I.%03d',stnum,sernum,subim);
  a = readim(str,[64 64],header);
  tc(im) = sum(sum(a(xloc:xloc+xsz,yloc:yloc+ysz)));
  ser(im) = sernum;
end
pe(slnum) = std(tc)/mean(tc)*100;
pe2(slnum) = std(detrend(tc))/mean(tc)*100;
mn(slnum) = mean(tc);
end
subplot(211)
plot(mn)
subplot(212)
plot(pe2)

return
str = sprintf('study %05d, series %03d, slice %d, err (%dx%d) = %.2f',stnum,stser,slnum,xsz,ysz,pe);
title(str);
ylabel('normalized image intensity')
xlabel('image number/time (s)');
