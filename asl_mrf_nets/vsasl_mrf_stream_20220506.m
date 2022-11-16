function vsasl_mrf_stream_20220506(inputFile, doRecon, netDir)
% function vsasl_mrf_stream(inputFile, doRecon, netDir)
% reconstruction, realignment and segmentation of tissue types
% then calculation of parameter estimates
% using Neural Network Regression

if doRecon
    fprintf('start recon ...')
    !rm *.nii
    sprec1_3d(inputFile,'fy','l','N', 'n', 128);
    % figure out the name of the NIFTI file with the data
    vname = dir('vol*.nii');
    figure
    lightbox(vname(1).name);
    [data h] = read_nii_img(vname(1).name);
    
    
    % realignment of time series
    !rm realign*
    fprintf('\nRealigning with MCFLIRT ... ')
    str = ['!mcflirt -in ' vname(1).name ' -cost normmi -out realigned -smooth 3 -bins 128 -refvol ' num2str(h.dim(5)-1) ]
    eval(str)
    !gunzip realigned.nii.gz
    [data h] = read_nii_img('realigned.nii');
    
    figure
    plot(mean(data,2))
    title('Time course of Spatial Average')
    drawnow
    
    % Do the clustering algorithm to segment the images
    Ntypes = 3;
    Niter = 50;
    fprintf('\nClustering segmentation with %d tissue Types and %d iterations ... ', Ntypes, Niter);
    timecourses = mrf_segment(vname(1).name, Ntypes, Niter);
    save segment_timecourses.txt timecourses -ascii
    drawnow

else
        [data h] = read_nii_img(inputFile);
end
%%
%
fprintf('\nCalculating perfusion with Neural Net ... ');
% calculate perfusion from time series using NN regression
% Network, weights and final scaling factor are stored in a mat file
str = [netDir '/Mynet_cbf.mat']
load(str)

cbf = NNreg_imgseries_20220506(data , Mynet, scale, 0.5);
cbf = cbf*6000; % adjust units to ml/min/100g

% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting perfusion file ... ');
write_nii('CBF_NNRE.nii', cbf(:) , h,0);
figure
subplot(211)
cbf = lightbox('CBF_NNRE.nii', [0 100],[]);
subplot(212)
hist(cbf(:), 100)
drawnow
%%
fprintf('\nCalculating T1 with Neural Net ... ');
str = [netDir '/Mynet_r1.mat']
load(str)

r1 = NNreg_imgseries_20220506(data , Mynet, scale, 0.5);
t1 = 1./r1;
t1(isinf(t1)) = 0;
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting T1 file ... ');
write_nii('T1_NNRE.nii', t1(:) , h,0)
figure
subplot(211)
lightbox('T1_NNRE.nii', [0 3],[]);
subplot(212)
hist(t1(:), 100);
drawnow

%%
fprintf('\nCalculating T2 with Neural Net ... ');
str = [netDir '/Mynet_r2.mat']
load(str)

r2 = NNreg_imgseries_20220506(data , Mynet, scale, 0.5);
t2 = 1./r2;
t2(isinf(t2)) = 0;
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting T2 file ... ');
write_nii('T2_NNRE.nii', t1(:) , h,0)
figure
subplot(211)
lightbox('T2_NNRE.nii', [0 3],[]);
subplot(212)
hist(t2(:), 100);
drawnow

%%
fprintf('\nCalculating CBV with Neural Net ... ');
str = [netDir '/Mynet_cbv.mat']
load(str)

cbv = NNreg_imgseries_20220506(data , Mynet, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting CBV file ... ');
write_nii('CBV_NNRE.nii', cbv(:) , h,0);
figure
subplot(211)
cbv=lightbox('CBV_NNRE.nii', [0 0.02],[]);
subplot(212)
hist(cbv(:), 100)
drawnow

%%
fprintf('\nCalculating BAT with Neural Net ... ');
str = [netDir '/Mynet_bat.mat']
load(str)

bat = NNreg_imgseries_20220506(data , Mynet, scale, 0.5);
% fix the header so that it writes out a single frame of floats
h.dim(5)  = 1;
h.datatype = 16;
h.bitpix = 32;

fprintf('\nWriting BAT file ... ');
write_nii('BAT_NNRE.nii', bat(:) , h,0);
figure
subplot(211)
bat = lightbox('BAT_NNRE.nii', [0 0.5],[]);
subplot(212)
hist(bat(:), 100)
drawnow
%%
!mkdir mrf_results
!rm Con* *map* log* spmJr*
!mv segment*.txt  *.nii *.img *.hdr mrf_results
