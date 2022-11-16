spnum = 7;
nsl = 12;
slnum = 7;
offset = nsl + slnum + 1;
header = 0;
numim = 600;
xloc = 32;
yloc = 32;
sz = 6;
sz2 = sz-1;

im = 1;
str = sprintf('sp%d/sl%d.phs.%03d',spnum,slnum,im);
a = readim(str,[64 64],header);
a1 = ones(size(a));
a1(xloc:xloc+sz2,yloc:yloc+sz2) = 1.2.*ones([sz sz]);
show(a.*a1);
str = sprintf('spiral %d, slice %d',spnum,slnum);
title(str);
pause
for im = 1:numim;
  str = sprintf('sp%d/sl%d.phs.%03d',spnum,slnum,im);
  a = readim(str,[64 64],header);
  tc(im) = sum(sum(a(xloc:xloc+sz2,yloc:yloc+sz2)));
end
freq = tc./(sz*sz*1000*pi*0.025);
plot(freq)
str = sprintf('UM3T(CERD-U/D) spiral run %d, slice %d, err = %.2f',spnum,slnum,pe);
title(str);
ylabel('frequency (arb offset)')
xlabel('image number/time (s)');
