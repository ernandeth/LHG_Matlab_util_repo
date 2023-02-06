function im = cartesian_recon3d(cart_data, dim, smaps)
%
% function im = cartesian_recon3d(cart_data, dim [,smaps])
%
% very simple 3D FFT of cartesian k-space data.
% combines coils using the sum of squares
% 
%     cartdata: is assumed to be in a 2D matrix that has dimensions 
%             Nk x Ncoils
%     dim : vector specifying the 3D dimensions of the k-space data
%

im=zeros(dim);
Ncoils = size(cart_data,2);
alldims = [dim Ncoils];
if nargin==3
    fprintf('\nCartesian Recon using using sensitivity maps\n')
    % scaling factor:  sum of the sensitivity maps
    SM = sum( abs(smaps).^2, 4) ;
    SM(isinf(SM)) = 0;
    SM(isnan(SM)) = 0;
    
% code from recon3dflex
%     im = div0( sum( conj(args.smap) .* im, 4), ...
%            sum( abs(args.smap).^2, 4) );
%     im = im .* (sum(abs(args.smap),4) > args.tol);
        

else
    smaps = ones(alldims);
    SM = 1;
end

for c=1:size(cart_data,2)
    tmp = reshape(cart_data(:,c), dim);
    tmp(isnan(tmp)) = 0;

    tmp = fftshift(fftn(tmp));
    tmp = permute(tmp,[2 1 3]);

    % add the squares of the coil images
    if nargin< 3
        im = im + tmp.^2;
    else
        tmp = tmp .* conj(smaps(:,:,:,c));
        im = im +tmp;
    end

end

if nargin< 3
    im = abs(sqrt(im));
end
% remove voxels where sensitivity is too low
im(SM<1e-3) = 0;
im = im ./ SM;
        
%figure
% ov([],im, round(dim(1)/2), round(dim(2)/2),round(dim(3)/2),0); 
% lightbox(abs(im));
return
