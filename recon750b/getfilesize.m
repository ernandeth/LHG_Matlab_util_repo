function fileSize = getfilesize(strFile,fmt)
% EXAMPLE
% strFile = '/subjects/111130eb/111129_0826_trigECG.dat';
% fileSizeMB = getfilesize(strFile,'MB')
% fileSizeBytes = getfilesize(strFile); % default

s = dir(strFile);
sizeBytes = s.bytes;

if ~exist('fmt','var')
    fmt = 'bytes';
end

switch lower(fmt)
    case 'mb'
        fileSize = sizeBytes*9.536*1e-7;
    case 'bytes'
        fileSize = sizeBytes;
    otherwise
        error('fmrilab:invalidfmt','Invalid fmt "%s", currently only works for bytes or MB',fmt);
end