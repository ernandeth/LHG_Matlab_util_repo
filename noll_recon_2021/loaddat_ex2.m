function dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
% order of data in file (rev/for, views(spirals), slices, coils)

% $Id: loaddat_ex2.m 1266 2014-03-20 18:12:40Z klitinas $
% $HeadURL: svn+ssh://klitinas@anger.engin.umich.edu/svn/matlab/recon/trunk/loaddat_ex2.m $

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
nodd = rem((scaninfo.nphmult*scaninfo.npr),2);
%nodd = 0;
% point to right view
% Doug changed the next line to make npr > 1 work
coilskip = scaninfo.nslices*(concat+1)*(scaninfo.nphmult*npr+1+nodd)*4*scaninfo.ndat; 
nproj = slnum2*(concat+1)*(scaninfo.nphmult*npr+1+nodd) + phnum*scaninfo.npr + 1;
echostep = (scaninfo.nphmult*scaninfo.npr+1+nodd);
for coilnum = 1:ncoils
  coiloff = coilskip*(coilnum-1);
  offs = scaninfo.headersize + nproj*ndat*4 + coiloff;
  fseek(fid,offs,'bof');
  for lp = 1:npr
    ktmp = fread(fid,2*ndat,'int16');
    dat(:,lp+(coilnum-1)*(concat+1)*npr) = ktmp(1:2:end) + i.*ktmp(2:2:end);
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
