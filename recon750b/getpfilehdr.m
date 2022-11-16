function hdr = getpfilehdr(strFile)
% EXAMPLE
% strFile = './P32768.7';
% hdr = getpfilehdr(strFile);

% Check for file
if ~exist(strFile,'file')
    error('fmrilab:invalidfile','File %s does not exist',strFile);
end

% Read header
fid = fopen(strFile,'r','l');
hdr = read_gehdr(fid);

% If header has character shift, read it using modified read_gehdr
if hasheadershift(hdr)
    offset = 5;
    while  hasheadershift(hdr) && offset > -6
        hdr = read_shifted_gehdr(fid,offset);
        offset = offset - 1;
    end
end

% Close file
fclose(fid);