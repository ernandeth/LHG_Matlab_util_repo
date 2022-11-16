function [checksum,numBytes] = checksum(strFile)
% EXAMPLE
% strFile = '/work/new_scanner/sample_run_set/sample_pfile_02.7';
% [checksum,numBytes] = checksum(strFile)
% 
% Author - Krisanne Litinas
% $Id: checksum.m 580 2013-05-20 16:16:55Z klitinas $

% Do checksum with unix syntax
strCmd = sprintf('cksum %s',strFile);
[status,strResult] = unix(strCmd);

if status
    error('fmrilab:invalidunixcall','Error when executing command "%s"',strCmd);
end

% Get rid of any random carriage returns if necessary
strResult = regexprep(strResult,'\r\n|\n|\r','');

strPattern = '(?<checksum>\d*) (?<bytes>\d*) (?<strFile>.*)';
sMatch = regexpi(strResult,strPattern,'names');
checksum = str2double(sMatch.checksum);
numBytes = str2double(sMatch.bytes);