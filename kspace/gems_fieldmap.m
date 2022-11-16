% script to make a field map
im1=fidFT('gems_20131031_01.fid/fid');
im2=fidFT('gems_20131031_02.fid/fid');
deltaTE = 0.001; % seconds

msk = ones(size(im1));
msk(abs(im1)<1.5*std(abs(im1(:)))) = 0;
lightbox(msk);

gamma = 267.5 * 1e2  ; % rad/s/Gauss

delt = angle(im1./im2) ; % radians

fm_gauss =   msk.*delt/(gamma * deltaTE) ;
figure;
lightbox(fm_gauss,[-0.1 0.1],5);

save fm_gauss.mat fm_gauss