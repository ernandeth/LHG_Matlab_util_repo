function timecourses = aal_timecourses(name);
% function timecourses = aal_timecourses(name);
%
% this function extracts 116 timecourses from the specified file using the AAL
% ROIs.  it returns a matrix where each column is a time course from one of
% those ROIs
%
% it's very important that the AAL image and your images are in the same
% space.
%
% What we do is use spm reslice and generate a new AAL template by
% resampling the old one to match the size of the data we're working with.
%
% for every project, you need to generate a new rROI_MNO_V4.nii and rmyAAL_ROI.mat file
% that match the space you're working with.
%

load ~/matlab/SPM5/aal/ROI_MNI_V4_List.mat

[raw h] = read_img(name);

%{
j=nifti('rROI_MNI_V4.nii');

myAAL_ROI = j.dat(:,:,:);
figure; lightbox(myAAL_ROI);

save rmyAAL_ROI myAAL_ROI
%}
load ~/matlab/SPM5/aal/rmyAAL_ROI.mat

timecourses = zeros(h.tdim, length(ROI));

for n=1:length(ROI)
    inds = find(myAAL_ROI==ROI(n).ID);
    tmp = raw(:,inds);
    timecourses(:,n) = mean(tmp,2);
end

return
