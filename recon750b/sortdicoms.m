function sortdicoms(strDir,strDirOut)
% EXAMPLE
% strDirFlat = '/work/new_scanner/dicom_stuff/epi/a';
% strDirOut = '/work/new_scanner/dicom_stuff/epi/a_sorted';
% sortdicoms(strDirFlat,strDirOut);

% If unspecified, work on pwd
if ~exist('strDir','var')
    strDir = pwd;
end

% If unspecified, define output directory
[strParentPath,strName] = fileparts(strDir);
if ~exist('strDirOut','var')
    strDirOut = fullfile(strParentPath,[strName '_sorted']);
end

% Create output directory
if ~exist(strDirOut,'dir')
    mkdir(strDirOut)
end

% Get list of dicoms to read
casFiles = getpatternfiles('.*',strDir,'cas');

% Have to sort the files based on number (eg make it so z100 doesn't
% directly follow z10
casFiles = sortcasfiles(casFiles,'z');
numFiles = length(casFiles);

% Loop through and sort based on series number and description
iThisSeries = 1;
for i = 1:numFiles
    strFile = casFiles{i};
    s = dicominfo(strFile);
    seriesNum = s.SeriesNumber;
    strSeriesDesc = s.SeriesDescription;
    strSeriesDesc = lower(strrep(strSeriesDesc,' ','_'));
    strDirNameSeries = sprintf('series_%d_%s',seriesNum,strSeriesDesc);
    
    if isfield(s,'NumberOfTemporalPositions')
        numTimePts = s.NumberOfTemporalPositions;
    else
        numTimePts = 1;
    end
    numSlices = s.ImagesInAcquisition;
    numFilesThisSeries = numTimePts * numSlices;
    
    % Increment instance of this series if not at last image
    if iThisSeries <= numFilesThisSeries
        if iThisSeries == 1 % Make new directory
            strDirSeries = fullfile(strDirOut,strDirNameSeries);
            if exist(strDirSeries,'dir')
                strDirSeries = strrep(strDirSeries,'series','2_series');
            end
                mkdir(strDirSeries);    
        end
    end
    
    if iThisSeries < numFilesThisSeries
        iThisSeries = iThisSeries + 1;
    elseif iThisSeries == numFilesThisSeries
        iThisSeries = 1;
    end
    

    %fprintf('\nCopying "%s" to %s',strFile,strDirSeries);
    copyfile(strFile,strDirSeries) 
end