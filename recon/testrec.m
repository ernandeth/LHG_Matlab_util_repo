[args,scaninfo,kinfo] = sprec1('P16384.7','cp','s',20);
%[args,scaninfo,kinfo] = sprec1('P16384.7');


if args.norec == 1
  return
end

fid = fopen(args.pfile,'r','l');

nim = args.nim;
iph = args.iph;
isl = args.isl;
zflg = args.zflg;
mapdda = scaninfo.mapdda;
unfact = args.unfact;
nphases = scaninfo.nphases;
nslices = scaninfo.nslices;

if (args.rtype == 'cp') | (args.rtype == 'ft')
% for gridding recon
% set up kaiser-bessel functions
scfactor = 0.04;
OVER = 2; % oversampling factor
NW = 1024; % gridding width
kb_gridl = 1.5; % gridding width (in orig coords)
kb_beta = 6.7; % Jackson's betas 1.5 = 6.7, 2.0 = 9.1, 2.5 = 11.5
kbin = [-NW/2:NW/2]/NW;
kbweights = besseli(0,kb_beta*sqrt(1.0-(2.0*kbin).^2));
kx = real(kinfo.k)*OVER;
ky = imag(kinfo.k)*OVER;
%  /* kaiser-bessel correction */
wcmax = sinh(sqrt(kb_beta^2))/sqrt(kb_beta^2);
wctmp = (2*pi*kb_gridl*[-nim/2:nim/2-1]/nim/OVER).^2 - kb_beta^2 + eps;
wc1d = real(sqrt(wctmp)./sin(sqrt(wctmp))*wcmax);
wc = wc1d'*wc1d;
[xx yy] = meshgrid([-nim/2:nim/2-1]);
wc = wc.*(abs(xx + i.*yy) < 0.55*nim)*scfactor;
subim = [-nim/2:nim/2-1] + nim*OVER/2 + 1;
% set up apodization
apod = fermi(nim*OVER,max(abs(kinfo.k))*OVER,-1);
if (args.rtype == 'cp') 
  % set up standard temporal interpolators
  tmax = max(kinfo.t);
  % for now, make the max segment length 4 ms long - make this user selectable
  L = ceil(tmax/0.004);
  tau = tmax/L;
  t2 = [0:L]*tau;
  % for now, make the temporal interpolators min-max interpolators assuming 
  % a generic triangular histogram with a FWHM of 100 Hz
  df = 100;
  [ll llt] = meshgrid([0:L]);
  ggt2 = sinc((ll-llt)*tau*df).^2;
  gb2 = sinc((([0:L]'*tau*ones(size(kinfo.t)) - ones([L+1 1])*kinfo.t)*df)).^2;
  beta = 0; %.00001;
  % a2 has the 1D temporal interpolators
  a2 = inv(ggt2 + beta*eye(L+1))*gb2;
  % tints has the gridded versions of a2
  % make time map
  datgr = reshape(kbgrid_mex(kx,ky,kinfo.t,ones(size(kx)),nim*OVER,kb_gridl*OVER,kbweights(:)),[nim nim]*OVER);
  tmap = real(datgr)./(imag(datgr)+eps);
  tints = zeros(nim*OVER,nim*OVER,L+1);
  for l = 0:L
    tmpint = interp1(kinfo.t,a2(l+1,:),tmap);
    tints(:,:,l+1) = tmpint.*apod;
  end
  tints(find(isnan(tints))) = 0;
end
end

% for iter recon
% ??

if (args.map == 1)
  fm = zeros(nim,nim,nslices);
  stsl = 0; ensl = nslices; 
  stph = 0; enph = 2; iph = 0; args.unfact = 1; args.unf = 0;
  for slnum = stsl:ensl-1
    for phnum = 0:1
      % do field map recon here
      dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
        for coilnum = 1:scaninfo.ncoils 
          % need to fix this stuff for coils/for/rev/npr
          datin = (dat(end:-1:1).').*kinfo.kdens;
          datgr = reshape(kbgrid_mex(kx,ky,real(datin),imag(datin),nim*OVER,kb_gridl*OVER,kbweights(:)),[nim nim]*OVER);
          imtmp = ift2(datgr.*apod);
        if (phnum == 0)
          imkb = wc.*imtmp(subim,subim);
        else 
          imkb = wc.*imtmp(subim,subim) .* conj(imkb);
        end 
      end % for coilnum
    end % for phnum
    % add fancy smoother
    fmtmp = angle(imkb);
    fmmag = sqrt(abs(imkb));
    fmsm = angle(conv2(fmmag.*exp(i*fmtmp),ones([7 7]),'same'));
    fm(:,:,slnum+1) = fmsm/(scaninfo.mapdel*1e-6)/2/pi; % in Hz
  end % for slnum
  save fmfile fm;
else
    if (args.rtype ~= 'ft')
      load fmfile;
    end
    if (iph == 0) 
      stph = ((1-zflg)*mapdda*unfact); 
      enph = nphases*unfact; 
    else
      stph = iph-1+mapdda; 
      enph = iph+mapdda; 
    end
    if (isl == 0) 
      stsl = 0; ensl = nslices; 
    else   
      stsl = isl-1; ensl = isl; 
    end
    for phnum = stph:enph-1
      phnumout = phnum+1-mapdda*unfact;
      for slnum = stsl:ensl-1
        for coilnum = 1:scaninfo.ncoils 
         dat = loaddat_ex2(phnum,slnum,scaninfo,fid);
         % for reverse spiral
         datin = (dat(end:-1:1).').*kinfo.kdens;
         datgr = reshape(kbgrid_mex(kx,ky,real(datin),imag(datin),nim*OVER,kb_gridl*OVER,kbweights(:)),[nim nim]*OVER);
         if (args.rtype == 'ft')
           % for standard Fourier (gridding) recon
           imtmp = ift2(datgr.*apod);
           imkb = wc.*imtmp(subim,subim);
         else
           % for time segmented recon
           imtmp = ift2(squeeze(datgr.*tints(:,:,1)));
           imkb = imtmp(subim,subim);
           for nt = 1:L
             imtmp = ift2(datgr.*tints(:,:,nt+1));
             imkb = imkb + imtmp(subim,subim).*exp(-i*2*pi*fm(:,:,slnum)*nt*tau); 
           end
           imkb = imkb.*wc;
         end
         writeim(slname(slnum+1,phnum),flipud(abs(imkb)'));
         % show(abs(imkb))
         % [kkx kky] = meshgrid([-32:31]);
         % datgr = griddata(real(kinfo.k),imag(kinfo.k),datin,kkx,kky);
         % datgr(isnan(datgr)) = 0;
         % im = ft2(datgr);
        end % for coilnum = 1:args.ncoils 
      end % for slnum = stsl:ensl-1
    end %for phnum = stph:enph-1
end % if (args.map == 1) 

