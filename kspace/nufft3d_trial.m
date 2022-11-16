function [img, Gob] = nufft3d_trial(kxyz, kdata, fov, dims)
omega = zeros(size(kxyz));
for d=1:3
        omega(:,d) = 2*pi * kxyz(:,d) * fov(d) / dims(d);
end

N = dims;
mask = true(N);
nufft_args = {N,6*ones(size(N)),2*N,N/2,'table',2^10,'minmax:kb'};
Gob = Gnufft(mask, [{omega},nufft_args(:)']);

% kdata = Gob*phantom3d(64);
img = qpwls_pcg1([],Gob, 1, kdata(:), 1, 'niter', 40);
figure, im(reshape(img,dims(1),dims(2),dims(3)))

return