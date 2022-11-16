function maps = svd_segment(filename, Ncomps)
%function maps = svd_segment(filename, Ncomps)
%

if nargin<2
    Ncomps = 15;
end

TopPercentSigma = 80;

fprintf('\nFinding the top %d components in the noisiest pixels \n(pixels with top %d percent of the variance)\n...',Ncomps, TopPercentSigma)
[raw hdr] = read_img(filename);
if isfield(hdr,'originator')
    hdr=nii2avw_hdr(hdr);
end

sigma = std(raw,1);
msk = ccvarmask(sigma, TopPercentSigma);

figure;
subplot(211)
lightbox(reshape(sigma,hdr.xdim, hdr.ydim, hdr.zdim));
title('std. dev. map')
subplot(212)
lightbox(reshape(msk,hdr.xdim, hdr.ydim, hdr.zdim));
title('mask for analysis');

mraw = raw .* repmat(msk, hdr.tdim,1);
mraw(isinf(mraw))=0;
mraw(isnan(mraw))=0;

[u, s,v]=svd(mraw',0);

tcomponents = v(:,1:Ncomps);

spmJr(filename, abs(tcomponents), eye(Ncomps));

for n=1:Ncomps
    figure
    subplot(211)
    plot(tcomponents(:,n)), title (sprintf('component %d',n));
    set(gcf,'Name', 'SVD identified time components')
    subplot(212)
    str = sprintf('Zmap_%04d.img', n);
    lightbox(str,[3 10],[]);
end
return

function msk = ccvarmask(varimage , th)
% makes a mask that preserves the data with the top TH percentage of the
% values

ordered = sort(varimage(:));
Nintgrl = cumsum(ordered)/sum(ordered(:)) * 100;
thval = find(Nintgrl>100-th);
%subplot(211), plot(ordered); title('Ordered Std. Deviations')
%subplot(212), plot(Nintgrl); title('Intrageted , Normalized Std. Deviations')

thval = ordered(thval(1))
msk = varimage;
msk(msk < thval) = 0;
msk(msk>=thval) = 1;
return
