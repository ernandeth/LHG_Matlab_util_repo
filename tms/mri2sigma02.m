[scalp hscalp] = read_nii_img('~hernan/matlab/SPM5/templates/T1.nii');
[csf hcsf] = read_nii_img('~hernan/matlab/SPM5/apriori/csf.nii');
[w hw] = read_nii_img('~hernan/matlab/SPM5/apriori/white.nii');
[g hg] = read_nii_img('~hernan/matlab/SPM5/apriori/grey.nii');

scalp(find(csf>0)) = 0;
scalp(find(w>0)) = 0;
scalp(find(g>0)) = 0;

white = w/256 * 1.2 ;

grey = g/256 * 0.3;

csf = csf/256 * 1.6;

bone = scalp/256 * 0.01;

sigma = white + grey + csf + bone;

write_nii('MNI_sigma.nii', sigma*100, hscalp, 0)

save MNI_sigma.mat sigma

