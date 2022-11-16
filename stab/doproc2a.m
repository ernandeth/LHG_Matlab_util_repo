spnum = 1;
nsl = 12;
slnum = 6;
offset = nsl + slnum + 1;
header = 0;
numim = 300;
imoff = 300;
xloc = 32;
yloc = 29;
sz = 6;
sz2 = round(sz/2);
sz3 = sz-sz2-1;

asz = 20;
asz2 = round(asz/2);
asz3 = asz-asz2-1;

im = 1;
str = sprintf('sp%d/sl%d.%03d',spnum,slnum,im);
a = readim(str,[64 64],header);
a1 = ones(size(a));
a1(xloc-sz2:xloc+sz3,yloc-sz2:yloc+sz3) = 1.2.*ones([sz sz]);
show(a.*a1);
str = sprintf('spiral %d, slice %d',spnum,slnum);
title(str);
pause
imev = zeros([64 64]);
imod = zeros([64 64]);
numev = 0;
numod = 0;
for im = 1:numim;
  str = sprintf('sp%d/sl%d.%03d',spnum,slnum,im+imoff);
  a = readim(str,[64 64],header);
  if rem(im,2) == 1
    imev = imev + a;
    numev = numev+1;
  else
    imod = imod + a;
    numod = numod+1;
  end
  tc(im) = sum(sum(a(xloc-sz2:xloc+sz3,yloc-sz2:yloc+sz3)));
end
avdif = (imev/numev - imod/numod);
avim = (imev/numev + imod/numod);
corav = avdif(xloc-asz2:xloc+asz3,yloc-asz2:yloc+asz3);
spatstd = std(corav(:))*sqrt((numev+numod)/(sz*sz))./mean(mean(avim(xloc-asz2:xloc+asz3,yloc-asz2:yloc+asz3)))*100

pe = std(tc)/mean(tc)*100
pedet = std(detrend(tc))/mean(tc)*100
plot(tc./mean(tc))
str = sprintf('UM3T(CERD-U/D) spiral run %d, slice %d, err = %.2f',spnum,slnum,pe);
title(str);
ylabel('normalized image intensity')
xlabel('image number/time (s)');
