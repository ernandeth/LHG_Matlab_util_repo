function [img, Gob] = nufft2d(kxy, kdata, fov, dims)
% function [img, Gob] = nufft2d(kxy, kdata, fov, dims)
%
% a wrapper for the conjugate phase NUFFT recon
%
% kxyz is the K space trajectory:  (N x 2)
% kdata is the K space data : (Nx 1)
% fov is the field of view's size (1 x 2 )
% dims is the output image dimensions: (1 x 2)
%
% returns the complex image and the Gmri object (Gob)
%

% prepare the IRT toolbox
if ~exist('irtdir')
    run ~/matlab/irt/setup.m
end

% shouchang's
% weighted least squares recon:
% omega is the scaled version of ks (trajectory) from -pi to pi


omega = zeros(size(kxy));
for d=1:2
        omega(:,d) = 2*pi * kxy(:,d) * fov(d) / dims(d);
end
% make weighting function proportional to distance to center of
% kspace
wi =  eps + (omega(:,1).^2 + omega(:,2).^2 ).^0.5;
%wi(:) = 1;

N = dims;
mask = true(N);
nufft_args = {N,...
    6*ones(size(N)),...
    2*N,...
    N/2,...
    'table',...
    2^10,...
    'minmax:kb'};

Gob = Gnufft(mask, [{omega}, nufft_args(:)']);
img = Gob' * (wi .* kdata(:));
%img = qpwls_pcg1([],Gob, 1,  kdata(:), 1, 'niter', 20);
img = reshape(img, dims);

% Note:
% Gob = Gmri(kxyz, ig.mask, 'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args), with the units of kxyz = 1/meters and the unit of FOV = meters
% is equivalent to 
% kxyz = kxyz*0.24, Gob = Gnufft(mask, [{omega},nufft_args(:)']); (the unit of kxyz is 1/FOV and FOV is not specified).

% alternatively:  regularized
% R = Reg1(ig.mask,'beta',2^2,'type_penal','mat','offsets','3d:26','pot_arg',{'huber',0.05}); %% 0.9
% img = pwls_pcg1([],Gob, 1, kdata(:), R, 'niter',20); % 1

% fix the mask, zero-fill , etc.
%img = embed(img, ig.mask);

return

