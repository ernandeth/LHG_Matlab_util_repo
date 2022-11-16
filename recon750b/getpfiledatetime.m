function [strScanDate,strScanTime] = getpfiledatetime(strFile)
% EXAMPLE
% strFile = '/work/subjects/111129eb/111129_P01536.7';
% [strDate,strTime] = getpfiledatetime(strFile);
% 
% $Id: getpfiledatetime.m 598 2013-05-20 16:16:55Z klitinas $

% Open file
fid = fopen(strFile, 'r', 'l');
fseek(fid,16,-1);

% Read date and time bytes
strScanDate = fread(fid, 10, '*char')';
strScanTime = fread(fid, 5,'*char')';

% Date tacks on an extra one
strScanDate(end-2) = '';
strScanDate(end) = '';

% Close file
fclose(fid);