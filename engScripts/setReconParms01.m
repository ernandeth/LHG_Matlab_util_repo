load ReconParms


args
kinfo
scaninfo
fprintf('\n These are the Imaging params for Spiral Recon.  Are they OK? ');
fprintf('\n Change what you want, and type "return" when you are done.');
fprintf('\n (you probably just need scaninfo.nslices, scaninfo.opfov ...)')
keyboard

%
%% generate k-space trajectories
%
fprintf('\n Generating the new k-space trajectories ....');
kfact1 = -exp(i*pi*scaninfo.rotation/2);
kfact2 = exp(i*pi*args.numrot/2)/args.zoomer;
if (scaninfo.gtype == 0)
    % standard spiral
    [g,k,t,s]=spiralgradlx3(scaninfo.opfov,scaninfo.opxres,scaninfo.gts,scaninfo.slewrate,scaninfo.fsgcm,scaninfo.npr);
    vd = ones(size(k));
else
    % var dens spiral
    [g,k,t,s,vd]=dogradvd(scaninfo.opfov,scaninfo.opxres,scaninfo.gts,scaninfo.slewrate,scaninfo.fsgcm,scaninfo.npr,scaninfo.densamp);
    %  [g,k,t,s,vd]=dogradvd(scaninfo.opfov,scaninfo.opxres,scaninfo.gts,scaninfo.slewrate,scaninfo.fsgcm,2,900);
end
if (scaninfo.transpose == 0)
    k = imag(k) + i*real(k);
else
    k = -k;
end
k = k.*kfact1;
% for advancing samples
if (args.samp1 < 0)
    k = [ones([1 -args.samp1])*k(1) k];
    vd = [ones([1 -args.samp1])*vd(1) vd];
end
lk = length(k);
if lk < scaninfo.ndat
    k = [k k(lk)*ones([1 (scaninfo.ndat-lk)])];
    vd = [vd vd(lk)*ones([1 (scaninfo.ndat-lk)])];
else
    k = k(1:scaninfo.ndat);
    vd = vd(1:scaninfo.ndat);
end
lk = length(k);
% density comp factor
g = [0 diff(k)];
kdens = -abs(g).*sin(angle(g)-angle(k)).*vd;

% data modulation for shifting and bulk off-resonance
t = [0:lk-1]*scaninfo.ts;
xsh = args.lrshift/args.nim + scaninfo.pixshifth;
ysh = args.tbshift/args.nim - scaninfo.pixshiftv;
kmod = exp(i*2*pi*(args.pt*t +  xsh*real(k) + ysh*imag(k)));
k = k*kfact2;

kinfo.k = k;
kinfo.kdens = kdens;
kinfo.kmod = kmod;
kinfo.t = t;
%%
fprintf('\n Updated K-space trajectories and saved as ReconParms.mat .');
fprintf('\n The previous ReconParms.mat has been saved to  ReconParms_old.mat.');

!mv ReconParms.mat ReconParms_old.mat
save ReconParms.mat args scaninfo kinfo