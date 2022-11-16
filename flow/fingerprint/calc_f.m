PLD=1.8;
LabelDur=1.8;
SI_PD=abs(pCASL_ims(:,:,7)).*(1/(1-exp(-4/1.4)));
mean_sub_im=mean(abs(abs(pCASL_ims(:,:,9:2:end))-abs(pCASL_ims(:,:,10:2:end))),3);
CBF=(6000.*0.9.*(mean_sub_im)*exp(PLD/1.65))./(2*0.85*1.65.*SI_PD.*(1-exp(-LabelDur/1.65)));
CBF(CBF>200)=0;

figure;imagesc(CBF);colormap hot;
mask=roipoly();