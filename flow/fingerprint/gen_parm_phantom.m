% Generate Paramater maps based on a template of white matter/grey matter
% segmentation.a

w=lightbox('/home/hernan/matlab/spm8/apriori/white.nii');
g=lightbox('/home/hernan/matlab/spm8/apriori/grey.nii');
c=lightbox('/home/hernan/matlab/spm8/apriori/csf.nii');



t1map = 900 * w + 1400*g + 2800*c;
flowmap = 20 * w + 60*g + 0*c;
bat1map = 2 * w + 1.2*g + 0*c;
cbvmap = 0.005 * w + 0.02*g + 0*c;
bat2map = 0.5*bat1map;
b1map  = (pi/3)*ones(size(w));
ftmap  = 1.2*w + 0.6*g + 0*c;

slice = 40

t1map = t1map(1:2:end,1:2:end, slice);
flowmap = flowmap(1:2:end,1:2:end, slice);
bat1map = bat1map(1:2:end,1:2:end, slice);
bat2map = bat2map(1:2:end,1:2:end, slice);
b1map = b1map(1:2:end,1:2:end, slice);
cbvmap = cbvmap(1:2:end,1:2:end, slice);
ftmap = ftmap(1:2:end,1:2:end, slice);

subplot(321); lightbox(t1map); title('T1')
subplot(322); lightbox(ftmap); title('Flowthru Time')
subplot(323); lightbox(bat1map);title('BAT1')
subplot(324); lightbox(bat2map);title('BAT2')
subplot(325); lightbox(flowmap);title('flow')
subplot(326); lightbox(cbvmap);title('CBV')


figure
subplot(321); hist(t1map(:),50); title('T1')
subplot(322); hist(ftmap(:),50); title('Flowthru Time')
subplot(323); hist(bat1map(:),50);title('BAT1')
subplot(324); hist(bat2map(:),50);title('BAT2')
subplot(325); hist(flowmap(:),50);title('flow')
subplot(326); hist(cbvmap(:),50);title('CBV')
