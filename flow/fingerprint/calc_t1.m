%% Luis's Data
t_tag=load('D:\MRF\ASL\Data\151103dm\SR_schedule\nominal\t_tags.txt');
t_adjust=load('D:\MRF\ASL\Data\151103dm\SR_schedule\nominal\t_adjusts.txt');
t_delay=0.08;
t_aq=0.03;

t_samp2=t_tag+t_adjust+t_delay;
t_samp2=t_samp2(3:end);
t_samp=t_samp2(1:5:end);

temp=load_nii('vol_e2450_11_03_115_0059.nii');
pCASL_ims2=squeeze(double(temp.img));
% pCASL_ims2=double(abs(pCASL_ims(:,:,3:end)));
% clear pCASL_ims;
pCASL_ims(:,:,1)=mean(pCASL_ims2(:,:,3:5),3);
pCASL_ims(:,:,2)=mean(pCASL_ims2(:,:,6:10),3);
pCASL_ims(:,:,3)=mean(pCASL_ims2(:,:,11:15),3);
pCASL_ims(:,:,4)=mean(pCASL_ims2(:,:,16:20),3);
pCASL_ims(:,:,5)=mean(pCASL_ims2(:,:,21:25),3);
pCASL_ims(:,:,6)=mean(pCASL_ims2(:,:,26:30),3);
pCASL_ims(:,:,7)=mean(pCASL_ims2(:,:,31:35),3);
pCASL_ims(:,:,8)=mean(pCASL_ims2(:,:,36:40),3);
pCASL_ims(:,:,9)=mean(pCASL_ims2(:,:,41:45),3);
pCASL_ims(:,:,10)=mean(pCASL_ims2(:,:,46:50),3);
pCASL_ims(:,:,11)=mean(pCASL_ims2(:,:,51:55),3);
pCASL_ims(:,:,12)=mean(pCASL_ims2(:,:,56:60),3);


%% Katie's data

% t_tag=load('D:\MRF\ASL\Data\NVNC102015\SR_20\label_dur.txt');
% t_adjust=load('D:\MRF\ASL\Data\NVNC102015\SR_20\PAD.txt');
% t_delay=0.08;
% t_aq=0.02452;
% t_samp=t_tag+t_adjust+t_delay;
% t_samp=t_samp(3:end);
% load pCASL_ims_324;
% pCASL_ims=double(abs(pCASL_ims(:,:,3:end)));

%% ROI
T1=1;
A=1;
B=1;

figure;imagesc(pCASL_ims(:,:,end));
mask=roipoly();
for n=1:size(pCASL_ims,3)
    imtemp=pCASL_ims(:,:,n);
    t1curve(n)=mean(mean(imtemp(mask)));
end
t1curven=double(t1curve./max(t1curve));
x = lsqcurvefit(@fun_sr,[T1 A B], t_samp, squeeze(t1curven).');
x(1)
fit_t1curve=fun_sr(x,t_samp);

figure;plot(fit_t1curve);hold on;plot(t1curven);hold off;

%% Image

T1=1;
A=1;
B=1;

figure;imagesc(pCASL_ims(:,:,end));
mask=roipoly();
% pCASL_ims=pCASL_ims.*repmat(mask,[1 1 18]);
pCASL_ims=pCASL_ims.*repmat(mask,[1 1 12]);

t1curve=reshape(pCASL_ims,[size(pCASL_ims,1)*size(pCASL_ims,2),size(pCASL_ims,3)]);
t1curven=t1curve./max(t1curve(:));

for n=1:size(t1curven,1)
    x = lsqcurvefit(@fun_sr,[T1 A B], t_samp, double(squeeze(t1curven(n,:))).');
    x2(n,:)=x;
    fit_t1curve(n,:)=fun_sr(x,t_samp);
end
t1map=reshape(squeeze(x2(:,1)),[64 64]);
figure;imagesc(abs(t1map.*mask));caxis([0 3]);colormap hot;


