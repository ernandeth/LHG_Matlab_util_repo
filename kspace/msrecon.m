
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

addpath('/net/aronofsky/home/towi/matlab/recon');

clear all

fov=20;
xres=128;
nleaves=4;
gamp=2.2;
gslew=180;

% get spiral k-space trajectory, grads at 4e-6, RO at 5e-6
%  the interpolation here is of course wrong, we want to interpolate
%  the gradient waveform and calculate the k-space from that
[gx,gy,kx,ky,sx,sy]=genspi2(fov,xres,100e-3,4e-6,nleaves,gamp,gslew);

% interpolate the kspace traj. to match the readout samples
% (grads res:  4us,   sample res. = 5us.)
kxro=interp1([0:4e-6:4e-6*length(kx)-4e-6],kx,[0:5e-6:4e-6*length(kx)-5e-6]);
kyro=interp1([0:4e-6:4e-6*length(ky)-4e-6],ky,[0:5e-6:4e-6*length(ky)-5e-6]);
f=readraw('/net/aronofsky/data/towi/3Tnoise1/Pfrrun1.7',1);

%fix the length of the trajectory vector
kxro=kxro(1:length(f));
kyro=kyro(1:length(f));

kx1=kxro.';
ky1=kyro.';

%rotate the trajectory for the different shots
tmp=rot2d([kxro.' kyro.'],90);
kx2=tmp(:,1);
ky2=tmp(:,2);

tmp=rot2d([kxro.' kyro.'],180);
kx3=tmp(:,1);
ky3=tmp(:,2);

tmp=rot2d([kxro.' kyro.'],270);
kx4=tmp(:,1);
ky4=tmp(:,2);

% get the weighting
ww1=weight_vor([kx1],[ky1],nleaves);
ww2=weight_vor([kx2],[ky2],nleaves);
ww3=weight_vor([kx3],[ky3],nleaves);
ww4=weight_vor([kx4],[ky4],nleaves);



% get the acquired data
for index=1:1440
    fprintf('\rreconning....%d of 1440 ',index);
    f1=readraw('/net/aronofsky/data/towi/3Tnoise1/Pfrrun1.7',index);
    rr = rem(index,4);
    if rr==0, rr=4;end
    cmdstr = sprintf('im=k2image([kx%d],[ky%d],[f1],ww%d,xres);',rr, rr,rr);
    eval(cmdstr);
    im=imflip(imflip(im,1),2);
    if index>999 
        str = sprintf('sl1.%04d',index);
    else 
        str = sprintf('sl1.%03d',index);
    end
    
    writeim(str,round(1000*abs(im)),'short');
    
end

% %% field map part -- courtesy of Brad Sutton
% mapdel=2.5e-3;
% f1=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',1);
% f2=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',2);
% f3=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',3);
% f4=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',4);
% im1=k2image([kx1;kx2;kx3;kx4],[ky1;ky2;ky3;ky4],[ws1*f1;ws2*f2;ws3*f3;ws4*f4],ww,xres);
% f1=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',17);
% f2=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',18);
% f3=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',19);
% f4=readraw('/net/aronofsky/data/towi/3Tnoise1/Pmap2.7',20);
% im2=k2image([kx1;kx2;kx3;kx4],[ky1;ky2;ky3;ky4],[ws1*f1;ws2*f2;ws3*f3;ws4*f4],ww,xres);
% fm=angle(im1./im2)/mapdel;
% fm=real(imsmooth(fm,2));
% tt=[0:5e-6:5e-6*(length(kx1)-1)].';
% imc=k2image_we([kx1;kx2;kx3;kx4],[ky1;ky2;ky3;ky4],[ws1*f1;ws2*f2;ws3*f3;ws4*f4],ww,xres,fm,10,[tt;tt;tt;tt]);



% % the images may need to be flipped along both directions
% % -- at least to match D's recon

%subplot(221)
show(round(1000*abs(im))')
%subplot(222)
%show(round(1000*abs(im2))')
%subplot(223)
%show(round(1000*angle(im))')
%subplot(224)
%show(round(1000*angle(im2))')


%alternative recon method using griddata:
    
    % grid the spiral data
    [kx,ky] = meshgrid([-matrix_size/2: matrix_size/2-1]* (1/fov) );
    %keyboard
    kspace=griddata(real(kxy), imag(kxy), signal, kx, ky); 
    kspace(find(isnan(kspace))) = 0;

    rspace = fft2(kspace);
    
    