GAMMA = 26752  % rad/s/Gauss
SUB = 8;   % this is the factor by which we subsample the original data

% the tissue 
h = read_hdr('t1spgr.hdr');
object = read_img2(h,'t1spgr.img');
fov = h.xdim*h.xsize ;  % cm
fovz = h.zsize*h.zdim;

object = object(1:SUB:end, 1:SUB:end, 1:SUB:end);
h.xdim = h.xdim/SUB;
h.ydim = h.ydim/SUB;
h.zdim = round(h.zdim/SUB);
h.xsize = h.xsize*SUB;
h.ysize = h.ysize*SUB;
h.zsize = h.zsize*SUB;

NPIX = h.xdim*h.ydim*h.zdim;
object = reshape(object, h.zdim*h.xdim*h.ydim,1);
[x y z]  = meshgrid(...
    -fov/2:h.xsize:fov/2-1 , ...
    -fov/2:h.ysize:fov/2-1 , ...
    -fovz/2:h.zsize:fovz/2-1);
loc = [ reshape(x,NPIX,1) reshape(y, NPIX,1) reshape(z,NPIX,1) ];

fprintf('\nDone creating object and its coordinates\n');

% sensitivity patterns:
DECAY = 0.1;

sens1 = ones(h.xdim, h.ydim,h.zdim);
sens2 = sens1;
sens3 = sens1;
sens4 = sens1;

DIM = h.xdim;
for dist=1:DIM
    sens1(dist,:,:) = exp(-DECAY * abs(dist-10));
    sens2(:,dist,:) = exp(-DECAY * abs(dist-10));
    sens3(DIM-dist+1,:,:) = exp(-DECAY * abs(dist-10));
    sens4(:,DIM-dist+1,:) = exp(-DECAY * abs(dist-10));
end
sens1 = reshape(sens1, h.zdim*h.xdim*h.ydim,1);
sens2 = reshape(sens2, h.zdim*h.xdim*h.ydim,1);
sens3 = reshape(sens3, h.zdim*h.xdim*h.ydim,1);
sens4 = reshape(sens4, h.zdim*h.xdim*h.ydim,1);
sens = sens1.*sens2.*sens3.*sens4;

% apply the sensitivy pattern to the object
m1 = object.*sens1;
m2 = object.*sens2;
m3 = object.*sens3;
m4 = object.*sens4;

fprintf('\nDone applying sensitivity\n');



% 1 - Make some spiral waverforms:

kmax = (h.xdim/2)/fov;   % in-plane resolution
kzmax = (h.zdim/2)/fovz;
kzstep = kzmax/SUB;

% spiral gradient waveforms
%dt=0.5e-5;   % sample period.
dt=2e-5;   % sample period.
samp_time = 0.02;
t=[0:dt:samp_time];
t =t.^(2/3);

spfreq = 2*pi*(h.xdim/2) / t(end);

% build a spiral in a single plane
k_spiral=  complex(t .* sin(spfreq*t) , t.*cos(spfreq*t));
k_spiral = k_spiral * kmax/(abs(max(k_spiral)));   % (scaling)

% build the stack of spirals:
k=[];
for kz = -kzmax: kzstep: kzmax
    k = [k;  real(k_spiral)'  imag(k_spiral)'  (ones(size(k_spiral))*kz)'];
end

fprintf('\nDone making trajectory\n');

%%%%%%%%%%%%%%%%%%%%%%
% % use carteisan k space trajectory instead for checking
% kstep = kmax/16;
% [kx ky kz] = meshgrid(-kmax:kstep: kmax-kstep , ...
%     -kmax:kstep: kmax-kstep ,...
%     -kzmax:kzstep: kzmax-kzstep);
% k = [ reshape(kx,32*32*20 ,1) ...
%         reshape(ky,32*32*20 ,1) ...
%         reshape(kz,32*32*20 ,1)];
%%%%%%%%%%%%%%%%%%%%%%

% 2 - Compute the signal equation

kdata1 = sim_signal( loc , k , m1, [fov fovz]);
fprintf('\ncoil 1 signal done\n');
kdata2 = sim_signal( loc , k , m2, [fov fovz]);
fprintf('\ncoil 2 signal done\n');
kdata3 = sim_signal( loc , k , m3, [fov fovz]);
fprintf('\ncoil 3 signal done\n');
kdata4 = sim_signal( loc , k , m4, [fov fovz]);
fprintf('\ncoil 4 signal done\n');


NK = size(kdata1,1)*size(kdata1,2)*size(kdata1,3);
fprintf('\nDone with signal generation\n');

% this is the gradient waveform to put into the scanner in order to get
% the desired trajectory:
Gxy_spiral = diff(k_spiral / dt) * (2*pi)/GAMMA;
max_grad = max(real(Gxy_spiral));
Gz = diff(kz / dt) * (2*pi)/GAMMA;


% 3-  Do the reconstruction with straight FT
im1 = fft3d(kdata1);
im2 = fft3d(kdata2);
im3 = fft3d(kdata3);
im4 = fft3d(kdata4);

fprintf('\nDone with standard FFT recon - no sensitivity info \n');

% 4 - Now do the SENSE recon:
% encoding matrix including the sensitivity into the encoding:
save smallRes.mat

E = [];
for r = 1:NPIX
    rvec = loc(r,:);
    fprintf('\rlocation:  %2.2f %2.2f %2.2f -> %2.2f percent   ',rvec(1), rvec(2), rvec(3), r/NPIX*100);
    E = [ E ; sens1(r)*exp(-i*2*pi*k*rvec')  sens2(r)*exp(-i*2*pi*k*rvec' ) sens3(r)*exp(-i*2*pi*k*rvec' ) sens4(r)*exp(-i*2*pi*k*rvec' ) ]; 
end


% The recon is done by solving :  kdata = E * object
% ie - estimate the object by  inv(E'*E)* E * kdata

kdata = [reshape(kdata1,NK,1);
       reshape(kdata2,NK,1);
       reshape(kdata3,NK,1);
       reshape(kdata4,NK,1)];
        
im_est = pinv(E'*E)*E*kdata;
im_est = reshape(im_est, h.xdim, h.ydim, h.zdim);

