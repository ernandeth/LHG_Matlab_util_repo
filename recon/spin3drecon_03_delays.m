function [imgs, dims] = spin3drecon_03(Pfile, rotfile, gradfile)
% function [imgs, dims] = spin3drecon_03(Pfile, rotfile, gradfile)
%

dt = 4e-6; %seconds
gamma = 267.5 * 1e2 / (2*pi) ; % Hz/Gauss
MYDIM = 64;
FOV = 24;
dims = [MYDIM, MYDIM, MYDIM];

[dat nslices npoints nfids] = read_raw_3d(Pfile,0);
nframes = nfids/nslices;

dat = reshape(dat', npoints, nfids);
g = load(gradfile);
rots = load(rotfile);

k = gamma*cumsum(g)*dt;

% adjusting for acquisition delays
% k = k(1 : end-26-10,:);
% dat = dat(27 : end-10, :);

% in fact, we don't need these many points
k = k(1:2:end-50,:);
dat = dat(1:2:end-50, :);
%%

Kmax = 1.2*max((k(:)));
[kx , ky, kz ]= meshgrid( ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM), ...
    linspace(-Kmax,Kmax,MYDIM));

ctr = 1;



% % we can probably do it with just half the data
% inloc = inloc(1:2:end, :);
% indata = indata(1:2:end);
% %%%%%%%%%%%

gkdata = zeros(MYDIM, MYDIM, MYDIM);

imgs = [];

% first frame only:
f = 1;
% grab all of the shots ("slices")  each frame
beg = (f-1)*nslices + 1;
fin =  f*nslices;
tmpk = dat(: , beg:fin);
tmpk = tmpk(:)';
tic


load nbrs_weights.mat

for delay = 10:20
    
    kdel = k(1:end-delay , :);
    
    kxs = []; kys=[]; kzs=[];
    fixrot = [0 1 0;-1 0 0;0 0 1];
    %fixrot = [1 0 0;0 1 0;0 0 1];
    for(n=1:3:length(rots))
        m = rots(n:n+2,:);
        kr = (fixrot*m)*kdel';
        kxs = [kxs kr(1,:)'];
        kys = [kys kr(2,:)'];
        kzs = [kzs kr(3,:)'];
        ctr = ctr +1;
    end
    kxs = kxs(:)';
    kys = kys(:)';
    kzs = kzs(:)';
    
    
    tmpk_del = tmpk(1+delay:end);
    
    [gkdata nb_save wt_save] = myinterp3d(tmpk, kxs, kys, kzs, kx, ky, kz, 2);
    
    [gkdata] = myinterp3d_apply(tmpk_del,  kx, ky, kz, nb_save, wt_save);
    
    imbuf = fftshift(ifftn(fftshift(gkdata)));
    %im2 = fftshift(ifftn(fftshift(kim2)));
    toc
    fprintf('\n delay  %d : completed \n', delay);
    
    tmpdata = log(abs(gkdata));
    tmpdata(isinf(tmpdata)) = 0;
    
    subplot(211)
    lightbox(tmpdata);
    subplot(212)
    lightbox(abs(imbuf));
    drawnow
    pause(0.1)
    
    %figure
    %lightbox(log(abs(resamps(:,:,:,2))), [],[]);
    %figure
    imgs = [imgs; imbuf(:)'];
end

return

[imgs, dims] = spin3drecon_03_delays('P09216.7', 'rotmats.txt', 'grad.txt');

colormap hot

for d=1:size(imgs,1)
    tmp = imgs(d,:);
    tmp = reshape(tmp, dims(1), dims(2), dims(3));
    lightbox(abs(tmp(:,:, 25:40) ));
    drawnow
    pause(0.2);
end



