function nii2analyze(casFilesNii)
% EXAMPLE
% casFilesNii = {'./arun_01.nii;./rtmp_variance.nii;./run_01.nii};
% nii2analyze(casFilesNii);

% Loop through and convert each file
for i = 1:length(casFilesNii)
    strFile = casFilesNii{i};
    strFileBse = strrep(strFile,'.nii','');
    nii2avw(strFileBse)
end