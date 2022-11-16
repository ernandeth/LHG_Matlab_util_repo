function imgPCAfilter(rootname, range)
% function imgPCAfilter(rootname, range)
% this function will load up a time series of images
% it will use the pixels
fprintf('\nReading data file ...');
h = read_hdr(rootname);
data = read_img(rootname);
raw = data;

fprintf('\nCalculating mean image from time series ...');
meanimg = mean(raw,1);

fprintf('\nMasking Voxels in range %d', range);

data(meanimg < range(1), :) = 0;
data(meanimg > range(2), :) = 0;
mask = ones(size(meanimg));
mask (meanimg < range(1)) = 0;
mask (meanimg > range(2)) = 0;

fprintf('\nComputing SVD ...')
[U,S,V] = svd(data);

fprintf('\nBuilding regressors from top 10 components...')
X = U(1:10,:); 

fprintf('\Calculating parameter estimates of the components...');
beta_hats = pinv(X) * data;
fprintf('\Removing components  ...');
data = raw - beta_hats*X;

write_img_series('residual.img', h, data);
h.tdim = 1;
h.datatype = 1;
h.bitsperpixel = 1;
write_img('mask.img',h, mask);
