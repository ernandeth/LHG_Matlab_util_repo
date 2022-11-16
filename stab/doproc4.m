stnum = 36;
nsl = 12;
stser = 4;
offset = nsl + slnum + 1;
header = 7904;
numim = 600;
xloc = 31;
yloc = 29;
xsz = 3;
ysz = 3;

aaa = [];
slnum = 1;
for slnum=0:nsl-1
im = 1;
imnum = im*nsl + slnum;
sernum = floor(imnum/999)*20 + stser;
subim = imnum - floor(imnum/999)*999 + 1;
str = sprintf('%05d/%03d/I.%03d',stnum,stser,subim);
a = readim(str,[64 64],header);
aaa = [aaa a];
end
