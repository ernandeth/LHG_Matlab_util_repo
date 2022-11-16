load data4Matt
load fmap4Matt

deltat = 400/125000;
nx = 64; ny = 64; nz = 64;
kspace = [kxs; kys; kzs].'; % 1/cm
ti = linspace(0, deltat*size(kspace,1), size(kspace,1));
fov = [3 3 3]; % cm
mask = true(64,64,64);
N = size(mask);
args = {N, [6 6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
beta = 2^35; % regularization parameter

% no field map
% A = Gmri(kspace, mask, 'fov', fov, 'nufft_args', args);

% with field map
fmap = fmap(1:2:end,1:2:end,1:2:end);
A = Gmri(kspace, mask, 'fov', fov, 'nufft_args', args, ...
    'zmap', 1i*fmap, 'ti', ti, 'L', 6);

% finite differencing matrix
R = Reg1(mask, 'beta', beta);
R = R.C;

% conjugate gradient solver
x = A'*col(kdata.');
x = qpwls_pcg1(x, A, 1, kdata(:), R, 'niter', 30);
figure(1); im(embed(x, mask));

% construct a test ball
q = -nx/2:nx/2-1; y = -ny/2:ny/2-1; z = -nz/2:nz/2-1;
[xx, yy, zz] = ndgrid(q, y, z); clear q y z;
testball = zeros(size(mask));
testball(sqrt(xx.^2 + yy.^2 + zz.^2) <= nx/4) = 1;
balldat = A*testball;
figure(2); im(reshape(balldat, size(kdata)));
figure(3); im(kdata);