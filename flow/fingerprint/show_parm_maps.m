
load flows.mat flowmap
load trans.mat transmap
load disps.mat dispmap
load R_flows.mat Rmap
load kfors.mat kformap
load flipang.mat flipmap
load t1map.mat t1map
load cbva.mat volmap
load trans2.mat transmap2
%%
figure
subplot(421)
lightbox(t1map, [0 3],1);
title('T1 map  (s)');

subplot(422)
lightbox((flipmap), [ ], 1);
title('Flip Angle (deg)')

subplot(423)
lightbox(flowmap * 6000, [],1);
title('Perfusion in ml/100g/min')

subplot(424)
lightbox(volmap, [0 0.02], 1);
title('CBVa fraction');

subplot(425)
lightbox(transmap2);
title('Tissue BAT (s)');

subplot(426)
lightbox(transmap,[],1);
title('Arterial BAT (s)')

subplot(427)
lightbox(kformap, [],1);
title('MTR (1/s)')


subplot(428)
lightbox(Rmap, [0.8 1],1);
title('Goodness of fit (R)')

set(gcf,'Name',pwd)
set(gcf,'Position',[0 20 800 800])

colormap pink

print -depsc all_fp_estimates
print -djpeg all_fp_estimates



