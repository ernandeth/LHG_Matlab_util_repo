function casFiles = getallrecursivefiles(strDir)
% getallrecursivefiles.m - finds all files in all levels in a given directory
%
% INPUTS
% strDir - string, top level directory to search
%
% OUTPUTS
% casFiles - cell array of filenames found in input directory (all sublevels)
%
% EXAMPLE
% strDir = '/work/_workcopy/matlab';
% casFiles = getallrecursivefiles(strDir);

% $Id: getallrecursivefiles.m 589 2013-05-20 16:16:55Z klitinas $

% Get the data for the current directory
dirData = dir(strDir);

% Find which things are directories
dirIndex = [dirData.isdir];

% Get a list of the files
casFiles = {dirData(~dirIndex).name}';
if ~isempty(casFiles)
    casFiles = cellfun(@(x) fullfile(strDir,x),casFiles,'UniformOutput',false);  % Prepend path to files       
end

% Get a list of the subdirectories that aren't '.' or '..'
subDirs = {dirData(dirIndex).name};
validIndex = ~ismember(subDirs,{'.','..'});

% Loop over valid subdirectories
for iDir = find(validIndex)
    nextDir = fullfile(strDir,subDirs{iDir});    % Get the subdirectory path
    casFiles = [casFiles; getallrecursivefiles(nextDir)];  % Recursively call getallrecursivefiles
end

end
