function [datAll,ktmp] = loaddat_rds_coil(coilnum,slnum,scaninfo,fid)
% order of data in file (rev/for, views(spirals), slices, coils)

% Output is [npts x nframes] for a given slice/coil

% $Id: loaddat_rds_coil.m 1701 2015-04-13 16:53:50Z klitinas $

% if (scaninfo.sliceorder == 0)
%   slnum2 = round((slnum+1.5)/2) + rem(slnum,2)*round((scaninfo.nslices+0.5)/2) - 1;
% else
%   slnum2 = slnum;
% end
ndat = scaninfo.ndat;
npr = scaninfo.npr;

datAll = zeros(scaninfo.ndat,scaninfo.nslices);
for phnum = 0:scaninfo.nphases-1
    phaseOff =  scaninfo.ncoils * scaninfo.nslices * scaninfo.ndat * scaninfo.npr *4 * (phnum);
    sliceOff = scaninfo.ncoils*scaninfo.npr * scaninfo.ndat * 4 * (slnum);
    coilOff = scaninfo.ndat * 4 * coilnum;
    fseek(fid,phaseOff+sliceOff+coilOff,'bof');
    ktmp = fread(fid,2*ndat,'int16');
    datAll(:,phnum+1) = ktmp(1:2:end) + i.*ktmp(2:2:end);   
end