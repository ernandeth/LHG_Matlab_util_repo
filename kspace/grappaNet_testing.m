% load a data set:
%{
[data scaninfo]= read_raw_3d('/export/data/djfrey/spiral3d/20220614/ph_newaqtest5/P35328.7',0);
kxyz = load('/export/data/djfrey/spiral3d/20220614/ph_newaqtest5/ktraj_cart.txt');
plot(kxyz)
kv = load('/export/data/djfrey/spiral3d/20220614/ph_newaqtest5/kviews.txt');
%}
%

fprintf('\nReading data ....');
[data scaninfo]= read_raw_3d('/export/data/djfrey/spiral3d/20220822/human_SERIOS_4shot_defaults/P92672.7',0);
signal = squeeze(data(1,:,:));
kxyz    = load('/export/data/djfrey/spiral3d/20220822/human_SERIOS_4shot_defaults/ktraj_cart.txt');
kv      = load('/export/data/djfrey/spiral3d/20220822/human_SERIOS_4shot_defaults/kviews.txt');
%}


ndat = size(data,2) / size(kv,1);
ncoils = size(data,3)
nframes = size(data,1)
nechoes = size(data,2)/ndat

fprintf('\nFiguring out kspace trajectory and rotations ....\n');

nleaves = length(unique(kv(:,1)));
nslices = length(unique(kv(:,2)));
kxyzR = zeros(ndat,3,nleaves,nslices);
for leafn = 1:nleaves
    for slicen = 1:nslices
        R = kv((leafn-1)*nslices+slicen,end-8:end);
        R = reshape(R,[3 3])';
        tmpR = (R*kxyz')';
        
%         plot3(tmpR(:,1), tmpR(:,2), tmpR(:,3) );
%         axis([-2 2 -2 2 -2 2])
%         drawnow
%         pause(0.5)
        
        kxyzR(:,:,leafn,slicen) = tmpR;
    end
end
kxyzR = [reshape(kxyzR(:,1,:,:),[],1),reshape(kxyzR(:,2,:,:),[],1),reshape(kxyzR(:,3,:,:),[],1)];

kx = kxyzR(:,1);
ky = kxyzR(:,2);
kz = kxyzR(:,3);
ks = [kx(:) ky(:) kz(:)];


% Use synethetic data for troubleshooting
fprintf('\nInstead of real data, Use a numerical phantom and 6 arbitrary coils');
[sim_data b1]=  sim_kspace3d_sense(kx, ky, kz);
signal = sim_data;

%
%plot(kxyzR)


%%
CalRadius = 1; %<---- testing the gridding fucntion : use all of kspace
CalRadius =0.5


dim3=[1 1 1]*65;

% GNet = makeGrappaNet_20220729(kx, ky, kz, signal, dim3, CalRadius);
% GNet = makeGrappaNet_20220809(kx, ky, kz, signal, dim3, CalRadius);

[allNets rhos] = makeGrappaNet_20220909(kx, ky, kz, signal, dim3, CalRadius);

% Try to understand PSF
%[cart_psf dens] = GrappaNet_interpolate(allNets, ones(size(signal)), ks, floor(dim3));
%psf = cartesian_recon3d(cart_psf, floor(dim3));


[cart_data dens] = GrappaNet_interpolate(allNets, signal, ks, floor(dim3));

im = cartesian_recon3d(cart_data, floor(dim3));
%im = cartesian_recon3d(cart_data, floor(dim3), dens);

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
