% temp=t_samp.^(BAT).*exp(-t_samp.*BAT);
%temp=A.*(t_samp-BAT).^B.*exp(-(t_samp-BAT)./C);temp(temp<0)=0;

%% Luis's Data

t_tag=load('D:\MRF\ASL\Data\151103dm\bat_schedule\nominal\t_tags.txt');
t_adjust=load('D:\MRF\ASL\Data\151103dm\bat_schedule\nominal\t_adjusts.txt');
t_delay=0.08;
t_aq=0.03;
t_samp=t_tag(3:2:end)+t_delay;

temp=load_nii('vol_e2450_11_03_115_0049.nii');
pCASL_ims=double(temp.img);
pCASL_sub=abs(abs(pCASL_ims(:,:,3:2:end))-abs(pCASL_ims(:,:,4:2:end)));

h=fspecial('average');
for n=1:24,pCASL_sub(:,:,n)=imfilter(squeeze(pCASL_sub(:,:,n)),h);end
%% Katie's Data
 

% t_tag=load('D:\MRF\ASL\Data\NVNC102015\BAT_64\label_dur.txt');
% t_delay=load('D:\MRF\ASL\Data\NVNC102015\BAT_64\PLD.txt');
% 
% t_samp=t_tag(3:2:end)+t_delay(3:2:end);
% 
% 
% load pCASL_ims_323;
% pCASL_sub=abs(abs(pCASL_ims(:,:,3:2:end))-abs(pCASL_ims(:,:,4:2:end)));
% 

%% ROI

%Initial Guess
BAT=0;
A=1;
B=1;
C=1;

figure;imagesc(pCASL_sub(:,:,end));
mask=roipoly();
for n=1:size(pCASL_sub,3)
    imtemp=pCASL_sub(:,:,n);
    batcurve(n)=mean(mean(imtemp(mask)));
end
batcurve_mc=double((batcurve-mean(batcurve(1:3)))./max(batcurve));

[x fit_batcurve]= fun_bat2([BAT A B C], t_samp, squeeze(batcurve_mc).');
x
% fit_batcurve=fun_bat(x,t_samp);

figure;plot(t_samp,fit_batcurve,'.');hold on;plot(t_samp,batcurve_mc);hold off;

%% Image

%Initial Guess
BAT=0;
A=1;
B=1;
C=1;

figure;imagesc(pCASL_sub(:,:,end));
mask=roipoly();
pCASL_sub=pCASL_sub.*repmat(mask,[1 1 24]);

batcurve=reshape(pCASL_sub,[size(pCASL_sub,1)*size(pCASL_sub,2),size(pCASL_sub,3)]);
batcurve_mc=(batcurve-repmat(mean(batcurve(:,1:3),2),[1 24]))./max(batcurve(:));

% figure;
for n=1:size(batcurve_mc,1)
    [x fit_batcurve]= fun_bat2([BAT A B C], t_samp, double(squeeze(batcurve_mc(n,:))).');
    x2(n,:)=x;
%     fit_batcurve=fun_bat(x,t_samp);
% plot(t_samp,fit_batcurve);hold on;plot(t_samp,batcurve_mc(n,:));hold off;
% pause(0.05)
end
batmap=reshape(squeeze(x2(:,1)),[64 64]);
figure;imagesc(rot90(rot90(rot90(batmap.*mask))));caxis([0 10]);colormap hot;



