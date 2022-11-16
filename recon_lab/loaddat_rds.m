function dat = loaddat_rds(phnum,slnum,scaninfo,fid)
% order of data in file (rev/for, views(spirals), slices, coils)

% $Id: loaddat_rds.m 1482 2014-08-06 15:42:33Z klitinas $

if (scaninfo.sliceorder == 0)
  slnum2 = round((slnum+1.5)/2) + rem(slnum,2)*round((scaninfo.nslices+0.5)/2) - 1;
else
  slnum2 = slnum;
end
ncoils = scaninfo.ncoils;
ndat = scaninfo.ndat;
npr = scaninfo.npr;
concat = scaninfo.concat;
dat = zeros([ndat (concat+1)*npr*ncoils]);
% add one if odd numbered views
% nodd = rem((scaninfo.nphmult*scaninfo.npr),2);
%nodd = 0;
nodd = 0;
% point to right view

% Seek to start of this time frame
% phaseOff = scaninfo.ncoils * scaninfo.nslices * scaninfo.ndat * scaninfo.npr *4 * (phnum - 1);
phaseOff = scaninfo.ncoils * scaninfo.nslices * scaninfo.ndat * scaninfo.npr *4 * (phnum);
fseek(fid,phaseOff,'bof');

% Seek to slice
% sliceOff = scaninfo.npr * scaninfo.ndat * 4 * (slnum - 1);
sliceOff = scaninfo.ncoils*scaninfo.npr * scaninfo.ndat * 4 * (slnum);
fseek(fid,sliceOff,'cof');

% coilskip = scaninfo.nslices*(concat+1)*(scaninfo.nphmult+1+nodd)*npr*4*scaninfo.ndat; 
% nproj = slnum2*(concat+1)*(scaninfo.nphmult+1+nodd)*scaninfo.npr + phnum*scaninfo.npr + 1;
% echostep = (scaninfo.nphmult+1+nodd)*scaninfo.npr;
for coilnum = 1:ncoils
 
  for lp = 1:npr
    ktmp = fread(fid,2*ndat,'int16');
    dat(:,coilnum) = ktmp(1:2:end) + i.*ktmp(2:2:end);
  end
  if (concat == 1)
    offs = scaninfo.headersize + (nproj+echostep)*ndat*4 + coiloff;
    fseek(fid,offs,'bof');
    for lp = 1:npr
      ktmp = fread(fid,2*ndat,'int16');
      dat(:,lp+(coilnum-1)*(concat+1)*npr+npr) = ktmp(1:2:end) + i.*ktmp(2:2:end);
    end
  end
end
