%function coreg_series(workFile)
% Coregister a time series in Nifti format using the spm_coreg algorithm
% this allows for different cost functions and works better for 
% timeseries where the contrast is evolving 
%(eg - fingerprinting or ASL) with time
%
!rm vol*.mat rvol*
clear V
V = spm_vol(workFile);
refVol = V(1);
refMtx = V(1).mat;  % this will be the same for all the volumes (if it hasn't been realigned before)
tdim = size(V,1);
mat = zeros(4,4,tdim);
mat(:,:,1) = refMtx;
flags = [];
flags.cost_fun='nmi';
%flags.sep = [3 1];

rp = zeros(tdim,6);
oldmoveparms = zeros(1,6);
for m=[2:tdim]
    
    sourceVol = V(m);
    
    refVol = V(1);

    
    moveparms = spm_coreg(refVol, sourceVol, flags);
    
    m
    moveparms
    
    rp(m,:) = moveparms;
    
    M = inv(spm_matrix(moveparms))*refMtx
    mat(:,:,m) = M;
    V(m).mat = M;
end

%spm_write_4dvol(V, 'realigned_tmp.nii');

reslice_opts.which = 2;
reslice_opts.prefix='r';

[p nm e] = fileparts(workFile);
str = ['save ' nm '.mat mat']
eval(str)
str = ['save rp_' nm '.txt rp -ascii']
eval(str)

spm_reslice(workFile, reslice_opts);

return