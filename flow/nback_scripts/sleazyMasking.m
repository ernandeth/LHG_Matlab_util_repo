% make a mask based on the activation in the pre-training session:
[pre h] = lightbox('Zmap_0004.img');

[zth znth] = getFDRthreshold('Zmap_0004.img', 0.05);

msk =pre;
msk(abs(pre)<zth) =0;
msk(abs(pre)>=zth) =1;
lightbox(msk);
write_img('pretest_mask.img', msk, h);

z2=lightbox('Zmap_0002');
z2mask = z2.*msk;
lightbox(z2mask);
write_img('Zmap_0002masked.img', z2mask, h);

z=lightbox('Zmap_0001');
zmask = z.*msk;
lightbox(zmask);
write_img('Zmap_0001masked.img', zmask, h);
getFDRthreshold('Zmap_0001masked.img', 0.05);



% make a mask based on the activation in the post-training session:
[pre h] = lightbox('Zmap_0005.img');

[zth znth] = getFDRthreshold('Zmap_0005.img', 0.05);

msk =pre;
msk(abs(pre)<zth) =0;
msk(abs(pre)>=zth) =1;
lightbox(msk);
write_img('posttest_mask.img', msk, h);

z2=lightbox('Zmap_0002');
z2mask = z2.*msk;
lightbox(z2mask);
write_img('Zmap_0002masked2.img', z2mask, h);
getFDRthreshold('Zmap_0002masked2.img', 0.1);

z=lightbox('Zmap_0001');
zmask = z.*msk;
lightbox(zmask);
write_img('Zmap_0001masked2.img', zmask, h);
getFDRthreshold('Zmap_0001masked2.img', 0.1);