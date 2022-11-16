function modDate = getfiledatemodified(strFile,strFmt)
% EXAMPLE
% strFile = '/work/new_scanner/sample_pfile/sample_pfile.7';
% dateModified = getfiledatemodified(strFile)

strDir = fileparts(strFile);
sDir = dir(strDir);