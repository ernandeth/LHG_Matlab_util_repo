function [img, Gob] = nufft3d(kxyz, kdata, fov, dims, wi)
% function [img, Gob] = nufft3d(kxyz, kdata, fov, dims, weights)
%
% a wrapper for the conjugate phase NUFFT recon
%
% kxyz is the K space trajectory:  (N x 3)
% kdata is the K space data : (Nx 1)
% fov is the field of view's size (1 x 3 )
% dims is the output image dimensions: (1 x 3)
% weights; density compensation
%
% returns the complex image and the Gmri object (Gob)
%

%run  ~/matlab/irt/setup.m

% 
% ig = image_geom('nx',dims(1), ...
%     'ny', dims(2), ...
%     'nz',dims(3) ,...
%     'fov', fov, ...
%     'offsets', 'dsp');



% conj. phase reconstruction done here:
%{
Gob = Gmri(kxyz, ig.mask, ...
    'fov', fov, ...
    'basis', {'rect'}, ...
    'nufft', nufft_args);
nufft_args = {dims, 6*ones(size(dims)), 2*dims, dims/2, 'table', 2^10, 'minmax:kb'};
wi_basis = wi ./ Gob.arg.basis.transform; % trick! undo basis effect
minmax(wi_basis)
img = Gob' * (wi_basis .* kdata(:));
%}

% shouchang's
% weighted least squares recon:
% omega is the scaled version of ks (trajectory) from -pi to pi
omega = zeros(size(kxyz));
for d=1:3
        omega(:,d) = 2*pi * kxyz(:,d) * fov(d) / dims(d);
end

% make weighting function proportional to distance to center of
% kspace
% wi = sqrt(omega(:,1).^2 + omega(:,2).^2 + omega(:,3).^2 );
% wi = wi .^ 0.8;
% wi = 0.0001 + wi/max(wi);

% Doug's density comp in sprec:
% kdens = -abs(g).*sin(angle(g)-angle(k)).*vd;
% usees velocity of kspace and the distance to the center.

% 1.20.21
% Make a density compensation function
%{
Rmin = sqrt(sum( (omega(20,:) - omega(1,:) ).^2 ));
dcf = zeros(size(kdata(:)));
for n=1:length(dcf)
    parfor m=1:length(dcf)
        dist(m) = sqrt(sum( (omega(n,:) - omega(m,:)).^2 ) );
    end
    dcf(n) = length( dist < Rmin);
end
wi = dcf.^2;  % quadratic?  
%}

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

%img = qpwls_pcg1([],Gob, 1,  kdata(:), 0, 'niter', 20);
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

