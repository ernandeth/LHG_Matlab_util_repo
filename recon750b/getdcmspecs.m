function [xPix,yPix,sliceThickness,rows,cols] = getdcmspecs(strFile,strDirOut)
% getdcmspecs.m - kind of kludgy fix for anatomy stream in case dcmtk
% doesn't work.
% 
% EXAMPLE:
% strFile = '/Users/klitinas/data/mr750/dicom/s3323/ginx/e12345s2i1';
% [xPix,yPix,sliceThickness,rows,cols] = getdcmspecs(strFile);

% Author - Krisanne Litinas
% $Id: getdcmspecs.m 279 2012-10-10 18:42:25Z klitinas $

% Check input
if ~exist(strFile,'file')
    error('fmrilab:buildanatomy:invaliddcmfile','File "%s" does not exist',strFile);
end

% Read dicom and get fields of interest
s = dicominfo(strFile);
pixSpacing = s.PixelSpacing;
xPix = pixSpacing(1);
yPix = pixSpacing(2);
sliceThickness = s.SliceThickness;
rows = s.Rows;
cols = s.Columns;

[~,strBase] = fileparts(strFile);
strLogFile = fullfile(strDirOut,[strBase '.tmp']);
fid = fopen(strLogFile,'w');
fprintf(fid,'pixx:%g',xPix);
fprintf(fid,'\npixy:%g',yPix);
fprintf(fid,'\nslicethickness:%g',sliceThickness);
fprintf(fid,'\nrows:%g',rows);
fprintf(fid,'\ncols:%g',cols);


fclose(fid);