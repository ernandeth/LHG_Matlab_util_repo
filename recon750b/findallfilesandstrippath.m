function [casFiles,casFilesFound] = findallfilesandstrippath(strDir)
% findallfilesandstrippath.m - get list of files in directory including extension
% 
% INPUTS
% strDir - string, directory to search
% 
% OUTPUTS
% casFiles - cell array of strings containing filenames
% 
% EXAMPLE
% strDir = '/export2/data/klitinas/subjects/111129eb/anatomy';
% casFiles = findallfilesandstrippath(strDir)

% Author - Krisanne Litinas
% $Id: findallfilesandstrippath.m 585 2013-05-20 16:16:55Z klitinas $

strPatAll = fullfile(strDir,'*');
casFilesFound = simpleglob(strPatAll);
casFilesFoundLowerSort = lower(sort(casFilesFound));
[~,casNames,casExt] = cellfun(@fileparts,casFilesFoundLowerSort,'uni',false);
casFiles = strcat(casNames,casExt);