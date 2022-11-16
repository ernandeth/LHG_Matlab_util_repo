function [tSNR sSNR] = ASL_snr(mask_threshold)
% function [tSNR sSNR] = ASL_snr(mask_threshold)
[mc h] = read_img('mean_con');
th = mask_threshold*std(mc(:));

msk = zeros(size(mc));
msk(find (mc>th)) = 1;

[ts h] = read_img('sub.img');
ms = mean(ts,1);  % mean signal inside brain
tsd = std(ts,[],1);  % temporal std dev map
sNoise = std( ms(:).* ~msk(:));  % std. dev over the outside of the brain

tsd = reshape(tsd,h.xdim, h.ydim, h.zdim);
ms = reshape(ms,h.xdim, h.ydim, h.zdim);
msk = reshape(msk,h.xdim, h.ydim, h.zdim);

subplot(223); lightbox(ms ./ tsd+eps, [-2 2],4); title('tSNR per pixel');

ms = ms.* msk;
tsd = tsd.* ~msk;

subplot(221); lightbox(ms); title('temporal mean (inside brain)');
subplot(222); lightbox(tsd); title('temporal std. dev. (outside brain)');

colormap bone

signal = mean(ms(:));
tNoise = mean(tsd(:));

tSNR = signal/tNoise
sSNR = signal/sNoise
set(gcf,'Name',pwd)
subplot(224)
axis off
text(0,1,['tSNR= ' num2str(tSNR)]);
text(0,0.5,['sSNR= ' num2str(sSNR)]);
drawnow
return
