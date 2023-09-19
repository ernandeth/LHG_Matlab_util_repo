% mrf LSQ results summary : make figures from the results
dirs = dir('*SOS_FSE');
orthov=0


for dnum=1:length(dirs)

    cd(dirs(dnum).name)
    
%%
%

   % Do the analysis:
   [r1 r2  b1err cbf cbv bat res] = img_mrf_lsq('im_mag.nii'); 
   save estimates.mat
    
%% 
%}   
    load estimates.mat
    
    figure(dnum)
    set(gcf,'Name',dirs(dnum).name , 'Position',[1 1 1500 800])
    
    subplot(171)
    tmp = lbview(1./r1(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 3])
    axis image ; axis xy ; axis off
    title('T1')
    
    subplot(172)
    tmp = lbview(1./r2(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 0.2])
    axis image ; axis xy ; axis off
    title('T2')
    
    subplot(173)
    tmp =lbview(b1err(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([-0.15 0.15])
    axis image ; axis xy ; axis off
    title('B_1 Err')
    
    subplot(174)
    tmp = lbview(cbf(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 0.015])
    axis image ; axis xy ; axis off
    title('CBF')
    
    subplot(175)
    tmp = lbview(cbv(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 0.03])
    axis image ; axis xy ; axis off
    title('CBV')
    
    
    subplot(176)
    tmp = lbview(bat(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 0.3])
    axis image ; axis xy ; axis off
    title('BAT')
    
    subplot(177)
    tmp = lbview(res(:,:,5:2:end-5),'nrows',10);
    imagesc(tmp)
    caxis([0 1e-2] )
    axis image ; axis xy ; axis off
    title('Resid.')
    
    colormap gray
    drawnow
    print -dpng estimates_summary
    
    %  now orthogonal views
    
    figure(dnum + 10)
    set(gcf,'Name',dirs(dnum).name , 'Position',[1 1 500 1200])
    
    subplot(711)
    orthoview(1./r1);
    caxis([0 3])
    axis image ; axis xy ; axis off
    title('T1')
    
    subplot(712)
    orthoview(1./r2);
    caxis([0 0.2])
    axis image ; axis xy ; axis off
    title('T2')
    
    subplot(713)
    orthoview(b1err);
    caxis([-0.15 0.15])
    axis image ; axis xy ; axis off
    title('B_1 Err')
    
    subplot(714)
    orthoview(cbf);
    caxis([0 0.015])
    axis image ; axis xy ; axis off
    title('CBF')
    
    subplot(715)
    orthoview(cbv);
    caxis([0 0.03])
    axis image ; axis xy ; axis off
    title('CBV')
    
    
    subplot(716)
    orthoview(bat);
    caxis([0 0.3])
    axis image ; axis xy ; axis off
    title('BAT')
    
    subplot(717)
    orthoview(res);
    caxis([0 1e-2] )
    axis image ; axis xy ; axis off
    title('Resid.')
    
    colormap gray
    drawnow
    print -dpng estimates_summary_ortho
    %%
    cd ..
    
    
end
