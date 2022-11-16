spnum = 8;
nsl = 12;
slnum = 7;
offset = nsl + slnum + 1;
header = 0;
numim = 600;
xloc = 32;
yloc = 29;
sz = 6;
sz2 = sz-1;

im = 1;
str = sprintf('sp%d/sl%d.%03d',spnum,slnum,im);
a = readim(str,[64 64],header);
a1 = ones(size(a));
a1(xloc:xloc+sz2,yloc:yloc+sz2) = 1.2.*ones([sz sz]);
show(a.*a1);
str = sprintf('spiral %d, slice %d',spnum,slnum);
title(str);
pause
for im = 1:numim;
  str = sprintf('sp%d/sl%d.%03d',spnum,slnum,im);
  b = readim(str,[64 64],header);
  tc(im) = sum(sum(b(xloc:xloc+sz2,yloc:yloc+sz2)));
end
pe = std(tc)/mean(tc)*100
pedet = std(detrend(tc))/mean(tc)*100
plot(tc./mean(tc))
str = sprintf('UM3T(CERD-U/D) spiral run %d, slice %d, err = %.2f',spnum,slnum,pe);
title(str);
ylabel('normalized image intensity')
xlabel('image number/time (s)');
