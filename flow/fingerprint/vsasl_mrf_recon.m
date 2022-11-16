%
fprintf('start recon ...')
!rm *.nii
sprec1_3d('P_test','fy','l','N');
%}

vname = dir('vol*.nii');
workFile =  vname(1).name;
%{
fprintf('\ndoing (SPM) realignment on ....%s\n', workFile);
% realignment with SPM
options.rtm=1;
spm_realign(workFile, options);
spm_reslice(workFile);
workFile = ['r' workFile ];
%}
figure
lightbox(workFile);
[data h] = read_nii_img(workFile);

%clip the last frame:
if size(data,1)==202
    data=data(1:end-1,:);
    h.dim(5) = 201;
    write_nii(workFile, data, h, 0);
end

figure
plot(mean(data,2)) 
title('Time course of Spatial Average') 
drawnow
%%
Ntypes = 3;
Niter = 50;
fprintf('\nClustering segmentation with %d tissue Types and %d iterations ... ', Ntypes, Niter);
mrf_segment(workFile, Ntypes, Niter)

%%
fprintf('\nCalculating perfusion with Neural Net ... ');
% calculate perfusion from time series using NN regression
% Network, weights and final scaling factor are stored in a mat file
load Mynet_cbf_8layer.mat
cbf=NNreg_imgseries(data , Mynet_cbf, scale, 0.5);

% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting perfusion file ... ');
write_nii('cbf_test.nii', cbf(:) , h,0);
figure
lightbox('cbf_test.nii', [0.005 0.02],[]);
drawnow
%%
fprintf('\nCalculating T1 with Neural Net ... ');

load Mynet_r1.mat
r1 = NNreg_imgseries(data , Mynet_r1, scale, 0.5);
t1 = 1./r1;
t1(isinf(t1)) = 0;
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting T1 file ... ');
write_nii('T1_test.nii', t1(:) , h,0);
figure
lightbox('T1_test.nii', [],[]);
drawnow

%%
fprintf('\nCalculating CBV with Neural Net ... ');
load Mynet_cbv.mat
cbv = NNreg_imgseries(data , Mynet_cbv, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting CBV file ... ');
write_nii('CBV_test.nii', cbv(:) , h,0);
figure
cbv=lightbox('CBV_test.nii', [0 0.02],[]);
drawnow

%%
fprintf('\nCalculating BAT with Neural Net ... ');

load Mynet_bat.mat
bat = NNreg_imgseries(data , Mynet_bat, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting BAT file ... ');
write_nii('BAT_test.nii', bat(:) , h,0);
figure
bat = lightbox('BAT_test.nii', [0 0.02],[]);
drawnow