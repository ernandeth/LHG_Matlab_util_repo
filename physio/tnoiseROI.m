function tnoiseROI(data, header, DesMat)

% function tnoiseROI(data, header, DesMat)
%
% (C) 2007 Brandon W. Sur
% University of Michigan
% Last edit: July 12, 2007
%
% INPUTS:
%   data: uncorrected image data  
%   header: header from read_img
%   DesMat: Design Matrix
%   TR: TR is TR
%   Nframes: Number of time points in the data
% 
% OUTPUT:
%   This program identifies voxels with high temporal standard deviation 
%   and generates a file called 'tnoiseROI.nii', which contains only those voxels.


[m n] = size(data);

NewData = data; % NewData is used to store noise voxels

% 1.Calculate the variance of the original time series.

tSTD = zeros(n,1);

for k = 1:n
    tSTD(k,1) = std(data(1:m,k));
end

% 2.Based on the threshold chosen by user, construct the tSTD noise ROI. 

upper_frac = input('Enter the number for the fraction of voxels to be used for tnoiseROI(in percentage): ');

[sort_tSTD, I2] = sort(tSTD,'descend');

per = floor(upper_frac*n/100);
ind_noise_vox = I2(1:per,1);

[z3 z4] = size(ind_noise_vox);

for k = 1:n
    if k ~= ind_noise_vox(1:z3,1)
        NewData(:, k) = 0;
    end
end

% 3.Perform correlation calculation to remove any tnoise voxels with high correlation
% with design matrix.

[a, b] = size(DesMat); 

for k = 1:b
        [RHO PVAL] = corr(NewData(:,ind_noise_vox(k,1)), DesMat(:,k));
        if PVAL < 0.2 % exlcude voxels with a p-value less than 0.2
            NewData(:,ind_noise_vox(k,1)) = 0;
        end 
end

% 4.Generate a NIFTI file.
h = avw2nii_hdr(header);
write_nii('tnoiseROI.nii', NewData, h, 0);


