% sample recon scripts
 addpath ~/matlab/mfiles
 addpath ~/matlab/img
 addpath ~/rxdat/img
 addpath .
% P97792.7 - 4sh, ax A/P
% P98304.7 - 4sh, ax R/L
% P98816.7 - 4sh, sag S/I
% P99328.7 - 4sh, sag A/P
% P00000.7 - 4sh, cor S/I
% P00512.7 - 4sh, cor R/L
% P01536.7 - 12sh - sp-out
% P02048.7 - 12sh - sp-in

%sprec1('P02048.7','v','m')
%sprec1('P02048.7','v','t','2','h','C','0.5','d'','2')
sprec1('P61952.7','v','m','d','3')
sprec1('P61952.7','v','h','t','2','C','0.5','d','4','QI')
return

sprec1('P61952.7','v','h','t','2','C','0.5','d','-3','QI')
ainm3 = aread('vol_e476_10_19_121_0002');
sprec1('P61952.7','v','h','t','2','C','0.5','d','4','QI')
ain4 = aread('vol_e476_10_19_121_0002');
sprec1('P61952.7','v','h','t','2','C','0.5','d','5','QI')
ain5 = aread('vol_e476_10_19_121_0002');
sprec1('P61952.7','v','h','t','2','C','0.5','d','6','QI')
ain6 = aread('vol_e476_10_19_121_0002');

sprec1('P61952.7','v','h','0','C','0.5','d','2','QO')
sprec1('P61952.7','v','h','t','2','C','0.5','d','2','QO')
aout2 = aread('vol_e476_10_19_121_0002');
sprec1('P61952.7','v','h','t','2','C','0.5','d','4','QO')
aout4 = aread('vol_e476_10_19_121_0002');
sprec1('P61952.7','v','h','t','2','C','0.5','d','6','QO')
aout6 = aread('vol_e476_10_19_121_0002');

return
sprec1('P61952.7','v','m','QI')
sprec1('P61952.7','v','t','2','C','0.5','d','2')
%sprec1('P38400.7','v','t','1')
%sprec1('P62464.7','m','v','C','0.1')       % generate field map (ino fmfile.mat) 
%sprec1('P62464.7','h','v','C','0.1')       % generate field map (ino fmfile.mat) 
%sprec1('P38912.7','m','v','C','0.1')       % generate field map (ino fmfile.mat) 
%sprec1('P24576.7','d','-5','t','2','v','C','0.1')
%sprec1('P18944.7','m','v')
%sprec1('P18944.7','h','m','v')
%sprec1('P18944.7','h','v')
return

aa = zeros(size(xsall));
aa(:,:,:,1) = aread('vol_e9788_08_15_121_0001');
aa(:,:,:,2) = aread('vol_e9788_08_15_121_0002');
aa(:,:,:,3) = aread('vol_e9788_08_15_121_0003');
aa(:,:,:,4) = aread('vol_e9788_08_15_121_0004');
aa(:,:,:,5) = aread('vol_e9788_08_15_121_0005');


tsnrsp = mean(aa,4)./(std(aa,0,4)+eps);
tsnrgm = mean(xsall,4)./(std(xsall,0,4)+eps);

% sprec1('P38912.7','m','v','C','0.1')       % generate field map (ino fmfile.mat) 
% sprec1('P38912.7','h','t','2','v','C','0.1')       % generate field map (ino fmfile.mat) 
% figure; ashow('vol_e18_06_17_121_0002')
% title('150_m1')
% %sprec1('P40448.7','m','v','C','0.1')       % generate field map (ino fmfile.mat) 
% %sprec1('P40448.7','h','t','2','v','C','0.1')       % generate field map (ino fmfile.mat) 
% sprec1('P41472.7','m','v','C','0.1')       % generate field map (ino fmfile.mat) 
% %sprec1('P41472.7','h','t','2','v','C','0.1')       % generate field map (ino fmfile.mat) 
% sprec1('P38912.7','h','t','2','v','C','0.1')       % generate field map (ino fmfile.mat) 
% figure; ashow('vol_e18_06_17_121_0002')
% title('150_m2')
%sprec1('P62464.7','h','v','C','0.1')       % generate field map (ino fmfile.mat) 
%P38912.7 - 150
%P40448.7 - 180
%P41472.7 - 200
%P62464.7
%sprec1('P64000.7','v','h','m')    % redo field map using 1st field map to correct distortions
%sprec1('P64000.7','v','h','t','2','mat')  % recon time point 2 with field map corrections
%sprec1('P64000.7','v','h')    % recon all time points


