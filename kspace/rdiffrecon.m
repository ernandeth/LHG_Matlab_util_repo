
% Code Notes:
%  -  recon code provided by Sang-woo Lee and Brad Sutton
%  -  field map correction not implemented

% Pfile Notes:
% Pfile is organized as follows:
%  coil 1:
%  slice 1: baseline interleaves time1, interleaves time 2, etc...
%  slice 2: baseline interleaves time1, interleaves time 2, etc...
%  etc...
%  coil 2: same as above but some random skip between them
%  if #leaves and #points is odd then there is an extra baseline at end

% Notes on field map:
%   recon 2 different TE images (k2image)
%   divide them ./    (to avoid 2pi wraps)
%   take the angle and divide by delta_TE
%   may need to smooth but this is the field map
%   use k2image_we to get the phs corrected image
%    (use L=8 or 10, smoothing width 6)

% To recon CORRECTLY aliased images as for UNFOLD
% supply all the shots and use zeros in the data for
% the other interleaves

clear all

ppath='./';
%pfile='P02560.7';
nleaves=1;
imno=2;
slno=2;
doeachleaf=1;
gtype=1;
densamp=300;
sensefact=1.0;
fov=28;
xres=64;
xresout=64;
nims=32;
nsl=12;
domap=0;
swin=tri([-xresout/2:xresout/2-1],xresout/2,0); swin=swin'*swin;

% recon info
gamp=2.2;
gslew=180;
f0=readraw([ppath,pfile],1);
frsize=length(f0);
mappos=(nims*nleaves+1)*(slno-1);
impos=mappos+((imno-1)*nleaves+1);

% get spiral k-space trajectory, grads at 4e-6, RO at 5e-6
%  the interpolation here is of course wrong, we want to interpolate
%  the gradient waveform and calculate the k-space from that
if (gtype==1),
  [gx,gy,kx,ky]=genspivd(fov,xres,1000e-3,4e-6,densamp,sensefact*nleaves,gamp,gslew);
else,
  [gx,gy,kx,ky,sx,sy]=genspi2(fov,xres,1000e-3,4e-6,nleaves,gamp,gslew);
end;
kxro=interp1([0:4e-6:4e-6*length(kx)-4e-6],kx,[0:5e-6:4e-6*length(kx)-5e-6]);
kyro=interp1([0:4e-6:4e-6*length(ky)-4e-6],ky,[0:5e-6:4e-6*length(ky)-5e-6]);
kxro=kxro(1:frsize).';
kyro=kyro(1:frsize).';

% rotate leaves if necessary and get data
tt=[0:5e-6:5e-6*(frsize-1)].';
for m=1:nleaves,
  ws(m)=1;
  f(:,m)=readraw([ppath,pfile],impos+m-1);
  if (domap),
    fmap(:,m)=readraw([ppath,pfile],mappos+m-1);
  else,
    fmap(:,m)=zeros([frsize 1]);
  end;
  tmp=rot2d([kxro kyro],(m-1)*(360/nleaves));
  kxx(:,m)=tmp(:,1);
  kyy(:,m)=tmp(:,2);
end;
kxvec=kxx(:,1);
kyvec=kyy(:,1);
fvec=ws(1)*f(:,1);
fmapvec=fmap(:,1);
if (nleaves>1), 
  for m=2:nleaves,
    kxvec=[kxvec;kxx(:,m)];
    kyvec=[kyvec;kyy(:,m)];
    fvec=[fvec;ws(m)*f(:,m)];
    fmapvec=[fmapvec;fmap(:,m)];
    tt=[ tt ; [0:5e-6:5e-6*(frsize-1)].' ];
  end;
end;
kxxlrvec=kxx(1:densamp,1);
kyylrvec=kyy(1:densamp,1);
fflrvec=f(1:densamp,1);
if (nleaves>1), for m=2:nleaves,
  kxxlrvec=[kxxlrvec;kxx(1:densamp,m)];
  kyylrvec=[kyylrvec;kyy(1:densamp,m)];
  fflrvec=[fflrvec;f(1:densamp,m)];
end; end;

% get the weighting
ww=weight_vor(kxvec,kyvec,nleaves);
wwlr=weight_vor(kxxlrvec,kyylrvec,nleaves);

% now to recon
winwid=2.5;
im=k2image(kxvec,kyvec,fvec,ww,xresout,winwid);
if (doeachleaf),
  for m=1:nleaves,
    ffiltlr=zeros(size(kxxlrvec)); ffiltlr((m-1)*densamp+1:m*densamp)=ones([densamp 1]);
    ffilt=zeros(size(kxvec)); ffilt((m-1)*frsize+1:m*frsize)=ones([frsize 1]);
    ffilt2=zeros(size(kxvec)); ffilt2((m-1)*frsize+1:(m-1)*frsize+densamp)=ones([densamp 1]);
    imea(:,:,m)=k2image(kxxlrvec,kyylrvec,fflrvec.*ffiltlr,wwlr,xresout,winwid);
    imea2(:,:,m)=k2image(kxvec,kyvec,fvec.*ffilt2,ww,xresout,winwid);
    imeai(:,:,m)=k2image(kxvec,kyvec,fvec.*ffilt,ww,xresout,winwid);
    imea(:,:,m)=im_smooth(imea(:,:,m),0,1);
    imea2(:,:,m)=im_smooth(imea2(:,:,m),0,1);
    imeais(:,:,m)=im_filt(imeai(:,:,m),swin);
  end;  
end;
immap=k2image(kxvec,kyvec,fmapvec,ww,xresout,winwid);

%% field map part -- courtesy of Brad Sutton
if (domap),
  mapdel=2.5e-3;
  immap=k2image(kxvec,kyvec,fmapvec,ww,xresout,winwid);
  %fm=angle(im./immap)/mapdel;
  fm=angle(exp(-i*(angle(im)-angle(immap))))/mapdel;
  fm=real(imsmooth(fm,2));
  im2=k2image_we(kxvec,kyvec,fvec,ww,xresout,fm,winwid,tt);
  if (doeachleaf),
    % do not know how to implement a leaf field correction
  end;
else,
  im2=im;
end;

% other corrections here
im3=zeros(size(im));
im4=zeros(size(im));
im5=zeros(size(im));
im6=zeros(size(im));
im7=zeros(size(im));
mask=ones([xresout xresout]);
%mask=real((abs(imea(:,:,1))>(0.333*max(max(abs(imea(:,:,1)))))).*rect2d(96,96,16,16,128,128));
for m=1:nleaves,
  if (gtype==1),
    im3 = im3 + abs(imeai(:,:,m)).*exp(i*angle(imeai(:,:,m))).*exp(i*mask.*(angle(imea2(:,:,1))-angle(imea2(:,:,m))));
    im4 = im4 + abs(imeai(:,:,m)).*exp(i*angle(imeai(:,:,m))).*exp(i*(mean(mean(angle(imea(:,:,1))))-mean(mean(angle(imea(:,:,m))))));
    im5 = im5 + abs(imeai(:,:,m)).*exp(i*mask.*(angle(imeai(:,:,m))-angle(imea2(:,:,m))));
    im6 = im6 + abs(imeai(:,:,m)).*exp(i*mask.*(angle(imeai(:,:,m))-angle(imeais(:,:,m))));
    im7 = im7 + abs(imeai(:,:,m)).*exp(i*((-1)^(m+1))*mask.*(angle(imeai(:,:,m))-angle(imeais(:,:,m))));
  else,
    im3 = im3 + abs(imeai(:,:,m)).*exp(i*angle(imeai(:,:,m))).*exp(i*(mean(mean(angle(imeai(:,:,1))))-mean(mean(angle(imeai(:,:,m))))));
  end;
end;


%im=imflip(imflip(im,1),2);
%im2=imflip(imflip(im2,1),2);

%subplot(211)
show(round(1000*abs(im))')
%subplot(212)
%show(round(1000*abs(im2))')
%subplot(212)
%show(round(1000*angle(im))')
%show(round(1000*abs(im3))')
%subplot(224)
%show(round(1000*angle(im2))')

