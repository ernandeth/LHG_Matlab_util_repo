load 2018-09-27-GrantFigsDigiPathoPhantom.mat

close all

print -dpng cbva_pred
figure
lightbox(cbva_truth',r,[],[]); colormap copper
print -dpng cbva_tru
figure
plot(cbva_truth(:), cbva_pred(:),'*');
axis([r(1) r(2) r(1) r(2)]); axis square
xlabel('Truth'), ylabel('Estimate')
fatlines ; dofontsize(16)
print -dpng cbva_corr



r = [0.00 0.015]
lightbox(kfor_pred',r,[],[]);  colormap bone
print -dpng kfor_pred
figure
lightbox(kfor_truth',r,[],[]); colormap bone
print -dpng kfor_tru
figure
plot(kfor_truth(:), kfor_pred(:),'*');
axis([r(1) r(2) r(1) r(2)]); axis square
xlabel('Truth'), ylabel('Estimate')
fatlines ; dofontsize(16)
print -dpng kfor_corr


r = [0.40 1.3]
lightbox(r1_pred',r,[],[]);  colormap gray
print -dpng r1_pred

lightbox(r1_truth',r,[],[]); colormap gray
print -dpng r1_tru
figure
plot(r1_truth(:), r1_pred(:),'*');
axis([r(1) r(2) r(1) r(2)]); axis square
xlabel('Truth'), ylabel('Estimate')
fatlines ; dofontsize(16)
print -dpng r1_corr


r = [0 1.6]
lightbox(flip_pred',r,[],[]);  colormap jet
print -dpng flip_pred
figure
lightbox(flip_truth',r,[],[]); colormap jet
print -dpng flip_tru
figure
plot(flip_truth(:), flip_pred(:),'*');
r = [1.3 1.58]
axis([r(1) r(2) r(1) r(2)]); axis square
xlabel('Truth'), ylabel('Estimate')
fatlines ; dofontsize(16)
print -dpng flip_corr


%%
close all

for subject = 1:3
    
    switch subject
        case 1
            load 2018-05-09-GrantFigsISMRM.mat
            mkdir sub1
            cd sub1
        case 2
            load 2018-06-01-GrantFigsISMRM.mat
            mkdir sub2
            cd sub2
        case 3
            load 2018-06-03-GrantFigsISMRM.mat
            mkdir sub3
            cd sub3
    end
    figure
    r = [1 3]
    
    msk = bat_pred;
    msk(msk>0) = 1;
    se = strel('disk',1,4);
    msk = imerode(msk,se);
    lightbox(msk);
    
    lightbox((bat_pred .* msk)',r,[],[]);  colormap parula
    print -dpng bat_pred
    
    figure
    r = [0 120]
    lightbox((flow_pred.* msk)',r,[],[]);  colormap hot
    print -dpng flow_pred
    
    figure
    r = [0.00 0.02]
    lightbox((cbva_pred.* msk)',r,[],[]);  colormap copper
    print -dpng cbva_pred
    
    figure
    r = [0.00 0.025]
    lightbox((kfor_pred.* msk)',r,[],[]);  colormap bone
    print -dpng kfor_pred
    
    figure
    r = [0.20 1.2]
    lightbox(((1./r1_pred).* msk)',r,[],[]);  colormap gray
    print -dpng r1_pred
    
    figure
    r = [0 1.6]
    lightbox((flip_pred.* msk)',r,[],[]);  colormap jet
    print -dpng flip_pred
    
    cd ..
    
end