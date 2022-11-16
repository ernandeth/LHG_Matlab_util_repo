function casFiles = findallfiles(strDir)
% findallfiles.m - get list of files in directory including extension
% 
% INPUTS
% strDir - string, directory to search
% 
% OUTPUTS
% casFiles - cell array of strings containing filenames

% Author - Krisanne Litinas
% $Id: findallfiles.m 584 2013-05-20 16:16:55Z klitinas $

strPatAll = fullfile(strDir,'*');
casFilesFound = simpleglob(strPatAll);
casFilesFound = lower(sort(casFilesFound));
[~,casNames,casExt] = cellfun(@fileparts,casFilesFound,'uni',false);
casFiles = strcat(casNames,casExt);