function casFiles = getjustfiles(strDirParent)
% getjustfiles.m - returns list of files at top directory, skipping
% directories.
% 
% EXAMPLE
% strDirParent = '/export2/data/klitinas/subjects/111129eb';
% casFiles = getjustfiles(strDirParent);

% Author - Krisanne Litinas
% $Id: getjustfiles.m 595 2013-05-20 16:16:55Z klitinas $

s = dir(strDirParent);

numItems = length(s);
casFiles = {};
for i = 1:numItems
    si = s(i);
    blnDir = si.isdir;
    if ~blnDir && ~strcmpi(si.name(1),'.')
        casFiles = [casFiles;si.name];
    end
end