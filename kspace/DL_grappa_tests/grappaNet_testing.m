% edits% load a data set:
% editing 9.19.2023
% addpath /home/hernan/matlab/LHG_util_repo/kspace/DL_grappa_tests/


fprintf('\nReading data ....');
pfilename = dir('P*');
pfilename = pfilename(1).name

[raw,phdr] = readpfile(pfilename);
ndat = phdr.rdb.frame_size;
nechoes = phdr.rdb.user2;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
ntrains = phdr.rdb.user1;
nframes = phdr.rdb.user0;
tr = phdr.image.tr*1e-3;
dim = phdr.image.dim_X;
fov = phdr.image.dfov/10;

% reshape: ndat x (ntrains*nframes) x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils

raw = reshape(raw,ndat,ntrains,nframes,nechoes,ncoils);
% permute: ndat x ntrains x nframes x nechoes x ncoils
%           --> nframes x ndat x nechoes x ntrains x ncoils
raw = permute(raw,[3,1,4,2,5]);

% grab the first frame only for testing:
raw = raw(1,:,:,:,:);
nframes = 1;

% T2 compensation
%  We force the center of k-space to have the same magnitude in all echoes
for n=1:nframes
    for c=1:ncoils
        for s=1:nechoes
            for t=1:ntrains
                raw(n, :,s,t, c) = raw(n,:,s,t, c) / (abs(raw(n,end/2,s,t,c)));
                %plot(abs(signal2(:,s,c))); drawnow
            end
        end
    end
end
%
%

signal = reshape(raw, ndat*ntrains*nechoes , ncoils);
%
fprintf('\nFiguring out kspace trajectory and rotations ....\n');

% Process Trajectory (new)
% Load in kspace trajectory & view transformation matrices
ktraj = dir('ktraj*.txt');
ktraj = load(ktraj(1).name);
kviews = dir('kviews*.txt');
kviews = load(kviews(1).name);

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Allocate space for entire trajectory
ks = zeros(ndat,3,nechoes,ntrains);

% Transform each view
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;

        % Transform the trajectory
        ks(:,:,echon,trainn) = ktraj*T(:,:,mtxi)';
    end
end

ks = [reshape(ks(:,1,:,:),[],1),reshape(ks(:,2,:,:),[],1),reshape(ks(:,3,:,:),[],1)];


% Coil compression:  use SVD to reduce the number of effective coils
[comp_signal, S, CompressMat] = ir_mri_coil_compress(signal,'ncoil',10);

dim3=[1 1 1]*128;
CalRadius =0.5

%[allNets rhos] = makeGrappaNet_20230112(half_kx, half_ky, half_kz, half_signal, dim3, 0.9);
%[allNets rhos] = makeGrappaNet_20220909(half_kx, half_ky, half_kz, half_scaled_signal, dim3, CalRadius);
[allNets rhos] = makeGrappaNet_20230919(ks(:,1), ks(:,2), ks(:,3), comp_signal, dim3, CalRadius);

[cart_data dens] = GrappaNet_interpolate(allNets, comp_signal, ks, floor(dim3));

im_DLI = cartesian_recon3d(cart_data, floor(dim3));
im_DLI = im_DLI/norm(im_DLI(:));

orthoview(im_DLI);
title('DL interp. Full Data')

return
%%

% Use synethetic data for troubleshooting
%{
fprintf('\nInstead of real data, Use a numerical phantom and 6 arbitrary coils');
[sim_data b1]=  sim_kspace3d_sense(kx, ky, kz);

%}
%
%plot(kxyzR)
%%  Compare to model based with SENSE maps
%{
recon3dflex('frames',1, 'pfile',pfilename );
imc = readnii('coils_mag').*exp(1i*readnii('coils_ang'));
%smap = bart('ecalib -d0 -m1',fftc(imc,1:3));
tmp = zeros(size(imc));
for c=1:size(imc,4)
    tmp(:,:,:,c) = fft3d(squeeze(imc(:,:,:,c)));
end
smap = bart('ecalib -d0 -m1',tmp);

writenii('./smap_mag.nii',abs(smap));
writenii('./smap_ang.nii',angle(smap));

recon3dflex('smap',smap, 'frames',1 ,'pfile',pfilename );  
%recon3dflex( 'frames',1);  
!mv timeseries_mag.nii ref_timeseries_mag.nii

recon3dflex('smap',smap, 'frames', 1, 'clipechoes', [0 7], 'pfile',pfilename ); 
%recon3dflex('frames', 1, 'clipechoes', [0 7]); 
!mv timeseries_mag.nii half_ref_timeseries_mag.nii
%}


%%
% Remove second half of echo train
half_signal = signal(1:end/2,:);
half_kx = ks(1:end/2 , 1);
half_ky = ks(1:end/2 , 2);
half_kz = ks(1:end/2 , 3);
%%

%%
% repeat with sensitivity maps

%smap_2 = smap(:, :, end:-1:1,:);
smap_2 = smap(:,:,:,end:-1:1);
smap_2 = smap(:, end:-1:1,:,:);

for n=1:32,
    figure(1)
    lbview(1-(abs(smap_2(:,:,:,n))));
    im_DLI = cartesian_recon3d(cart_data(:,n), floor(dim3));
    figure(2)
    lbview(im_DLI);
    pause;
end

smap_2 = smap(:, end:-1:1,:,:);
cart_data_2 = cart_data;
%cart_data_2(:,23)= 10;
im_DLI_s = cartesian_recon3d(cart_data_2, floor(dim3), smap_2);

lbview(abs(im_DLI_s));

im_DLI_s = cartesian_recon3d(cart_data, floor(dim3), (smap));
im_DLI_s = im_DLI_s/norm(im_DLI_s(:));

title('DL interp. Full Data - SMAP')
writenii('DL_im_full', im_DLI);

%% Undersample the data R = 2
usignal = signal(1:end/2,:);
uks = ks(1:end/2,:);
[ucart_data dens] = GrappaNet_interpolate(allNets, usignal, uks, floor(dim3));

uim_DLI = cartesian_recon3d(ucart_data, floor(dim3));
uim_DLI = uim_DLI/norm(uim_DLI(:));

title('DL interp. Half Data')
writenii('DL_im_half', uim_DLI);

% repeat with sensitivity maps
uim_DLI_s = cartesian_recon3d(ucart_data, floor(dim3), smap);
uim_DLI_s = uim_DLI_s/norm(uim_DLI_s(:));
title('DL interp. Half Data -SMAP')



%% Undersample the data R = 4
u4signal = signal(1:end/4,:);
u4ks = ks(1:end/4,:);
[u4cart_data dens] = GrappaNet_interpolate(allNets, u4signal, u4ks, floor(dim3));

u4im_DLI = cartesian_recon3d(u4cart_data, floor(dim3));
u4im_DLI = u4im_DLI/norm(u4im_DLI(:));

% repeat with sensitivity maps
u4im_DLI_s = cartesian_recon3d(u4cart_data, floor(dim3), smap);
u4im_DLI_s = u4im_DLI_s/norm(u4im_DLI_s(:));

title('DL interp. Quarter Data')
writenii('DL_im_quart', u4im_DLI);

%% Try to understand PSF
[cart_psf dens] = GrappaNet_interpolate(allNets, ones(size(signal)), ks, floor(dim3));
psf = cartesian_recon3d(cart_psf, floor(dim3), smap);

%% Results
figure;
im_mirt = readnii('ref_timeseries_mag');
im_mirt = im_mirt/norm(im_mirt(:));
psf_mirt = readnii('psf.nii');
psf_mirt = psf_mirt/norm(psf_mirt(:));

orthoview(im_mirt);
title('MIRT. Full Data')

figure
uim_mirt = readnii('half_ref_timeseries_mag');
uim_mirt = uim_mirt/norm(uim_mirt(:));
orthoview(uim_mirt);
title('MIRT. Half Data')

% figure
% orthoview(psf_mirt);
% title('MIRT. PSF')

figure
orthoview(log(abs(psf)))
title('DLI PSF')
writenii('DL_psf', psf);

figure;
orthoview(im_DLI);
title('DLI Full Data')

figure
orthoview(uim_DLI);
title('DLI Half Data')

figure;
orthoview(im_DLI_s);
title('DLI (S-MAP) Full Data')

figure
orthoview(uim_DLI_s);
title('DLI (S-MAP) Half Data')

scl = [-1 1]*1e-2;

figure
orthoview(im_DLI-im_mirt)
title(sprintf('Difference - Full Data: RMS=%f', sqrt(sum((im_DLI(:)-im_mirt(:)).^2))));
caxis(scl);

figure
orthoview(uim_DLI-im_mirt)
title(sprintf('Difference - Half Data: RMS=%f', sqrt(sum((uim_DLI(:)-im_mirt(:)).^2))));
caxis(scl);

figure
orthoview(u4im_DLI-im_mirt)
title(sprintf('Difference - Quarter Data: RMS=%f', sqrt(sum((u4im_DLI(:)-im_mirt(:)).^2))));
caxis(scl);

figure
orthoview(im_DLI_s-im_mirt)
title(sprintf('Difference (S-MAP)- Full Data: RMS=%f', sqrt(sum((im_DLI_s(:)-im_mirt(:)).^2))));
caxis(scl);

figure
orthoview(uim_DLI_s-im_mirt)
title(sprintf('Difference (SMAP)- Half Data: RMS=%f', sqrt(sum((uim_DLI_s(:)-im_mirt(:)).^2))));
caxis(scl);

figure
orthoview(u4im_DLI_s-im_mirt)
title(sprintf('Difference - Quarter(SMAP) Data: RMS=%f', sqrt(sum((u4im_DLI_s(:)-im_mirt(:)).^2))));
caxis(scl);

mean(rhos)
median(rhos)
std(rhos)

%%
%
% code to eliminate redundant locations
%{
Ncenter=100;
Ndata = size(data,2);
Nechoes = size(kv,1);
inds = [];
repeats = [ndat/2 - Ncenter/2 : ndat/2 + Ncenter/2];

for echon = 1:Nechoes
    tmp = mean(kx(repeats));
    kx(repeats) = tmp;

    tmp = mean(ky(repeats));
    ky(repeats) = tmp;
    
    tmp = mean(kz(repeats));
    kx(repeats) = tmp;
    
    for f=1:nframes
        for c=1:ncoils
            tmp = mean(data(f,repeats, c),2 );
            data(f, repeats, c) = tmp;
        end
    end
    inds = [inds repeats(1:end/2-1) repeats(end/2+1:end)];
    repeats = repeats + ndat;
    
end
kx(inds) = [];
ky(inds) = [];
kz(inds) = [];
data(:,inds,:) = [];
%}
%{
fprintf('\nScaling all echoes to correct T2 effects...')
data2 = reshape(data, nframes, ndat, nechoes, ncoils);
for f=1:nframes
    for c=1:ncoils
        for e=1:nechoes
            data2(f,:,e,c) = abs(data2(f,end/2,1,c)) *data2(f,:,e,c) /  abs(data2(f,end/2,e,c));
        end
    end
end

data = reshape(data2,size(data));
%}
