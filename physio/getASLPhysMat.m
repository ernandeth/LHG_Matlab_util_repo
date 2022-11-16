function PM = getASLPhysMat(rootname, NumStd);
% function getASLPhysMat(rootname, NumStd);
% 
% identify CSF voxels as the brightest in the image (mean + Numstd * (standard
% deviations).   "mask,img" is written so you can see which pixels got
% chosen.
%
% extracts the mean time course in those voxels and builds a "Physio
% Matrix"
% 1- baseline regressor (column 1) 
% 2 - oscillating baseline for ASL (column 2)
% 3 - the mean centered, extracted time course (column 3)
% 4,5,6 - polynomial terms x, x^2, and x^3


[data, h] = read_img_series(rootname);
mdata = mean(data,1);
md = mean(mdata);
sd = std(mdata);
threshold = md+NumStd*sd;
inds = find(mdata> threshold);
ts = data(:,inds);
ts = mean(ts,2);
ts = ts - mean(ts);
ts = ts / sum(abs(ts));

x = [1:length(ts)] ;
x = x / sum(x) ;
x2 = x.^2 / sum(x.^2);
x3 = x.^3 / sum(x.^3);

baseline = ones(size(ts));
fbaseline = baseline;
fbaseline(1:2:end) = -1;
PM = [ baseline fbaseline ts x' x2' x3'];

mask = zeros(size(mdata));
mask(inds) = 1;

h.tdim = 1;
write_img('mask.img',mask,h);
write_hdr('mask.hdr',h);

return
