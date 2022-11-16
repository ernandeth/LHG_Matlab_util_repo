 
% SPM partitioned brain into these:

csf = read_img2('c3ht1spgr.nii');
 white = read_img2('c2ht1spgr.nii');
[ gr h] = read_img2('c1ht1spgr.nii');
 
% then we make a mask 
[whole nh] = read_nii_img('ht1spgr.nii'); 
tmp = whole;  tmp(whole>500)= 255; 
tmp(whole<=500)=0;

shell=tmp-white-gr-csf;
shell(shell<0) = 0;
write_nii('stuff.nii', shell,nh,0);
lightbox(shell);

% then we use FAST for the segmentation of the shell with 3 compartments
csf2 = read_img2('stuff_seg_0.nii');
bone = read_img2('stuff_seg_1.nii');
muscle = read_img2('stuff_seg_2.nii');

csf = csf + csf2 ;
lightbox(csf);

sigmaMap = (csf*1.73 + white*0.14 + gr*0.33 )/255 + bone*0.0042 + muscle* 0.33 ;
save sigmaMap sigmaMap
write_nii('sigmaMap.nii',sigmaMap*100,nh,0)