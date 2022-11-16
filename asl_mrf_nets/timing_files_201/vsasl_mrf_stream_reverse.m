function vsasl_mrf_stream(pfilename, doRecon)

if doRecon
    fprintf('start recon ...')
    !rm *.nii
    sprec1_3d(pfilename,'fy','l','N', 'n', 128);
end
%%
% figure out the name of the NIFTI file with the data
vname = dir('vol*.nii');
figure
lightbox(vname(1).name);
[data h] = read_nii_img(vname(1).name);

%fix error in P file: clip the last frame:
if size(data,1)==202
    data=data(1:end-1,:);
    h.dim(5) = 201;
    write_nii(vname(1).name, data, h, 0);
end

figure
plot(mean(data,2)) 
title('Time course of Spatial Average') 
drawnow
%
!rm realign*
fprintf('\nRealigning with MCFLIRT ... ')
str = ['!mcflirt -in ' vname(1).name ' -cost normmi -out realigned -smooth 3 -bins 128 -refvol ' num2str(h.dim(5)-1) ]
eval(str)
!gunzip realigned.nii.gz
[data h] = read_nii_img('realigned.nii');

%% Do the clustering algorithm to segment the images
Ntypes = 3;
Niter = 50;
fprintf('\nClustering segmentation with %d tissue Types and %d iterations ... ', Ntypes, Niter);
timecourses = mrf_segment('realigned.nii', Ntypes, Niter);
save segment_timecourses.txt timecourses -ascii
drawnow

%% clip the first 10 frames from the time series ?
% Load the shortened versions of the network
%data = data(10:end,:);

% 
fprintf('\nCalculating perfusion with Neural Net ... ');
% calculate perfusion from time series using NN regression
% Network, weights and final scaling factor are stored in a mat file
%load /home/hernan/matlab/asl_mrf_nets/Mynet_cbf_short.mat
%cbf = NNreg_imgseries(data , Mynet, scale, 0.5);

load /home/hernan/matlab/asl_mrf_nets/VSI_reverseOrder_211227/Mynet_cbf.mat
cbf = NNreg_imgseries(data , Mynet, scale, 0.5);

cbf = cbf*6000; % adjust units to ml/min/100g

% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting perfusion file ... ');
write_nii('CBF_test.nii', cbf(:) , h,0);
figure
subplot(211)
cbf = lightbox('CBF_test.nii', [0 100],[]);
subplot(212)
hist(cbf(:), 100)
drawnow

%%
fprintf('\nCalculating T1 with Neural Net ... ');

load /home/hernan/matlab/asl_mrf_nets/Mynet_r1_short.mat
r1 = NNreg_imgseries(data , Mynet, scale, 0.5);
t1 = 1./r1;
t1(isinf(t1)) = 0;
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting T1 file ... ');
write_nii('T1_test.nii', t1(:) , h,0)
figure
subplot(211)
lightbox('T1_test.nii', [0 3],[]);
subplot(212)
hist(t1(:), 100);
drawnow

%%
fprintf('\nCalculating CBV with Neural Net ... ');
load /home/hernan/matlab/asl_mrf_nets/Mynet_cbv_short.mat
cbv = NNreg_imgseries(data , Mynet, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting CBV file ... ');
write_nii('CBV_test.nii', cbv(:) , h,0);
figure
subplot(211)
cbv=lightbox('CBV_test.nii', [0 0.02],[]);
subplot(212)
hist(cbv(:), 100)
drawnow

%%
fprintf('\nCalculating BAT with Neural Net ... ');

load /home/hernan/matlab/asl_mrf_nets/Mynet_bat_short.mat
bat = NNreg_imgseries(data , Mynet, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting BAT file ... ');
write_nii('BAT_test.nii', bat(:) , h,0);
figure
subplot(211)
bat = lightbox('BAT_test.nii', [0 0.5],[]);
subplot(212)
hist(bat(:), 100)
drawnow
%%
!mkdir mrf_results
!rm Con* *map* log* spmJr*
!mv segment*.txt  *.nii *.img *.hdr mrf_results
