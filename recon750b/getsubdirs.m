function casDirs = getsubdirs(strDirParent)
% EXAMPLE
% strDirParent = '/export2/data/klitinas/subjects/111129eb';
% casDirs = getsubdirs(strDirParent);

% Author - Krisanne Litinas
% $Id: getsubdirs.m 600 2013-05-20 16:16:55Z klitinas $

s = dir(strDirParent);

numItems = length(s);
iDir = zeros(numItems,1);
% casItems = cell(numItems,1);
casDirs = {};
for i = 1:numItems
    si = s(i);
    blnDir = si.isdir;
    if blnDir && ~strcmpi(si.name(1),'.')
        casDirs = [casDirs;si.name];
    end
end

% casDirs = casItems(iDir > 0);