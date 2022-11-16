function dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
% order of data in file (rev/for, views(spirals), slices, coils)

if (scaninfo.sliceorder == 0)
  slnum2 = round((slnum+1.5)/2) + rem(slnum+1,2)*round((scaninfo.nslices+0.5)/2) - 1;
else
  slnum2 = slnum;
end

ncoils = scaninfo.ncoils;
ndat = scaninfo.ndat;
npr = scaninfo.npr;
concat = scaninfo.concat;

dat = zeros([ndat (concat+1)*npr*ncoils]);
% add one if odd numbered views
nodd = rem((scaninfo.nphmult*scaninfo.npr*(concat+1)),2);

% point to right view
nproj = slnum2*((concat+1)*scaninfo.nphmult*scaninfo.npr+1+nodd) + (concat+1)*phnum*scaninfo.npr + 1;

for coilnum = 1:ncoils
  coiloff = scaninfo.mcskip*(coilnum-1);
  offs = scaninfo.headersize + nproj*ndat*4 + coiloff;
  fseek(fid,offs,'bof');

  for lp = 1:((concat+1)*npr)
      % read a buffer and make it a complex array
      ktmp = fread(fid,2*ndat,'int16');
      dat(:,lp+(coilnum-1)*(concat+1)*npr) = ktmp(1:2:end) + i.*ktmp(2:2:end);
      
  end
end
