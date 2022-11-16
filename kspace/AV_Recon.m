%% NUFFT RECON Experiment
dt = 4e-6; %seconds
gamma = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma/2/pi ; % Hz/Gauss
% sl 0.8 for new sphere
%slthick = 0.6;
% myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/P08704.7';
% myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/P19456_big_crush';
%myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/spiral3d_G2.1_16slice_te055/P30208.7';
%myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/spiral3d_G2.1_16slice_te055_fov32/P31744.7';
%myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/spiral3d_G2.1_16slice_te055_fov32/P32768.7';
%myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/spiral3d_G2.1_16slice_te100/P31232.7';
myPfile = '/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/spiral3d_default_te054_nl4/P49664.7';

pfiles = dir('P*.7');
myPfile = pfiles(1).name;
% load the raw data:
% Nframes x Nslices*ndat x Ncoils
[raw scaninfo ] = read_raw_3d(myPfile ,0 );

nslices = scaninfo.nslices;
nechoes = scaninfo.npr;
nframes = scaninfo.nphases;
ndat = scaninfo.ndat;
ncoils = scaninfo.ncoils;
nfids = nechoes*nframes*nslices;
fov = scaninfo.opfov;
MYDIM = scaninfo.opxres;
slthick = scaninfo.slthick;
nleaves = scaninfo.npr;

plot(squeeze((abs(raw(1,:,:)))));
title('Echo train: Frame 1, all coils');
drawnow

fprintf('\nLoading grads...');
g = load('grad.txt');
glen = size(g,1);

% use the k-space trajectory calculated from the gradients,
% rather than the ideal trajectory
k2 = gamma*cumsum(g)*dt;
k2 = k2';
% must resample it to match the data acquisition rate
ks = zeros(3,ndat);
fprintf('\nResampling (down!) 2D trajectory to match signal acquisition rate...');
chunk = glen-ndat;
for n=1:3
    tmp = k2(n,:);
    %ftmp = fft(tmp);
    %ftmp = [ftmp(1:end/2-chunk/2) ftmp(end/2+chunk/2+1:end)];
    %ks(n,:) = real(ifft(ftmp)) * ndat/glen;
    ks(n,:) = interp1(linspace(0,1,glen), tmp, linspace(0,1,ndat));
end

Nramp = 500;
% throw out ramp
%ks = ks(:, Nramp+1:end-Nramp);
ks = ks';
% delays in gradient/acquisition?
Ndel = 55;
ks = [zeros(abs(Ndel),3); ks; zeros(abs(Ndel),3)];
ks(:,1) = circshift(ks(:,1),Ndel,1);
ks(:,2) = circshift(ks(:,2),Ndel,1);
ks = ks(abs(Ndel)+1:end-abs(Ndel),:);

ks=ks';

% if multiple interleaves, concatenate the rotated trajectories
if nleaves >1
    rotang = 2*pi/nleaves;
    ks_all = [];
    for n=0:nleaves-1
        c = cos(rotang*n); s=sin(rotang*n);
        R = [c s 0; -s c 0; 0 0 1];
        kstmp = R*ks;
        ks_all = [ks_all kstmp]; 

    end
    ks = ks_all;
end
    % throw out ramp
    ks = ks(:, Nramp+1:end-Nramp);

% filter to bring up dim stuff
% weights = abs(linspace(0,2,nslices))+1;
% 
% for s=1:nslices
%     ks(:,s) = ks(:,s) * weights(s);
% end

Kxymax = max(abs(ks(:)));
Kxymax = MYDIM/fov;  % Kmax
Kzmax = 1/(slthick)/2;
% Actually sampled  grid : stack of spirals
% kx and ky are always the same, but kz increases monotonically
% this was not recorded in the grad file
fprintf('\nCreating whole 3D sampled grid ...');
tmp = repmat(ks,  1, nslices);
kxs = tmp(1,:);
kys = tmp(2,:);
kzs = tmp(3,:);
ndat2 = ndat- Nramp*2;

kzz = linspace(-Kzmax, Kzmax,nslices);
for n=1:nslices
    %beg = ndat*(n-1)+1;
    beg = ndat2*(n-1)+1;
    %fin = ndat*n;
    fin = ndat2*n;
    kzs(beg:fin) = kzz(n);
end
 
ks = [kxs; kys; kzs];
ks=ks';
%kt = load('/Users/av/Documents/MATLAB/Spine_ASL/recon/spirals/ktraj_cart.txt');
raw2 = []; rawtmp = [];
for f = 1:nframes
    for c = 1:ncoils
        tmpk=raw(f,:,c);
        tmpk = reshape(tmpk,ndat*nleaves,nslices);
        % throw out ramp
        tmpk = tmpk(Nramp+1:end-Nramp,:);
        weights = abs(linspace(0,2,nslices))+1;
        for s=1:nslices
            tmpk(:,s) = tmpk(:,s) * weights(s);
        end
        tmpk(isnan(tmpk)) = 0;
        %rawtmp = [rawtmp, tmpk];
    end
    raw2(f, :) = reshape(tmpk, 1, []);
end

%%
ig = image_geom('nx', MYDIM, 'ny', MYDIM, 'nz',nslices ,'fov', [fov, fov, nslices*slthick], 'offsets', 'dsp', ...
    'mask', 'all-but-edge');
%ig = image_geom('nx', MYDIM, 'ny', MYDIM, 'nz',nslices ,'fov', [fov, fov, nslices*slthick], 'offsets', 'dsp');
%ig.mask = ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
N = ig.dim;

% Nspiral = 400;
% omega = linspace(0, 10*2*pi, Nspiral)';
% omega = pi*[cos(omega) sin(omega)] .* omega(:,[1 1])/max(omega);


% omega = atan2(ks(:, 2), ks(:, 1));
% omega = pi*[cos(omega) sin(omega)].*omega(:, [1 1])/max(omega);
omega = zeros(size(ks));
fovs = [ig.fov, ig.fov, ig.zfov];
for id=1:length(N)
	dx = fovs(id) / N(id);
	omega(:,id) = ks(:,id)*(2*pi)*dx;
end

omega = 2*pi*[ks(:,1) ks(:,2)] / MYDIM; % 
%omega = 2*pi*ks*[fov, fov, nslices*slthick]'./N;
%omega = mri_trajectory_stack(omega, nslices);
wi = abs(omega(:,1) + 1i * omega(:,2));
wi = wi / fov^2; % approximate scaling
wi = wi / (nslices*slthick); % Z dim scaling

nufft_args = {N, 6*ones(size(N)), 2*N, N/2, 'table', 2^10, 'minmax:kb'};

% Gn = Gnufft(ig.mask, {omega, nufft_args{:}});

Gm = Gmri(ks, ig.mask, ...
		'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);


%% Conjugate Phase Reconstruction
minmax(Gm.arg.basis.transform)
wi_basis = wi ./ Gm.arg.basis.transform; % trick! undo basis effect
minmax(wi_basis)
% 

 for n=2 % :nframes
        im_all = 0;
        for c = 1:4:ncoils
            tmpk = raw2(n,:,c)';
            xcp = Gm' * (wi_basis .* tmpk);
            xcp = embed(xcp, ig.mask);
            im(xcp, 'Conj. Phase Recon'), cbar
            im_all = abs(im_all) + i*abs(xcp);
        end
 end
im(im_all), cbar
%% Regularizer
% beta = 2^-21 * size(omega,1); % good for quadratic 'rect'
% R = Reg1(ig.mask, 'beta', beta);
% niter = 10;
% 
%     
% 
%  for n=2 % :nframes
%         im_all = 0;
%         for c = 1:ncoils
%             tmpk = raw2(n,:,c)';
%             xpcg = qpwls_pcg1(xcp(ig.mask), Gm, 1, tmpk, R.C, 'niter', niter);
%             xpcg = ig.embed(xpcg);
%             im(xpcg, '|x| pcg quad'), cbar
%             im_all = abs(im_all) + i*abs(xpcg);
%         end
%  end
% %im(im_all), cbar
% im(xpcg, ['|x| pcg quad, Delay: ', num2str(Ndel)]), cbar
%%
function omega = mri_trajectory_stack(omega2, N3)
o3 = single([0:(N3-1)]/N3 - 0.5)*2*pi;
o3 = repmat(o3, nrow(omega2), 1); % [N12,N3]
omega = repmat(omega2, N3, 1); % [N12*N3,2]
omega = [omega o3(:)]; % [N12*N3,3]
end
    
    
