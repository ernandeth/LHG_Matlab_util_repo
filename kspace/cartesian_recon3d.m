function im = cartesian_recon3d(cart_data, dim, dens)
%
% function im = cartesian_recon3d(cart_data, dim)
%
% very simple 3D FFT of cartesian k-space data.
% combines coils using the sum of squares
% 
%     cartdata: is assumed to be in a 2D matrix that has dimensions 
%             Nk x Ncoils
%     dim : vector specifying the 3D dimensions of the k-space data
%

im=zeros(dim);
if nargin==3
    fprintf('using density compensation')
    dens =reshape(dens, dim).^0.9;
else
    dens = ones(dim);
end

for c=1:size(cart_data,2)
    tmp = reshape(cart_data(:,c), dim) ./ dens;
    tmp(isnan(tmp)) = 0;
    
    figure
%    ov([],log(abs(tmp)), round(dim(1)/2), round(dim(2)/2),round(dim(3)/2),0); 
    subplot(211)
    lightbox(log(abs(tmp(:,:,1:3:end))));
%     
     tmp = fftshift(ifftn(tmp));
    subplot(212)
    lightbox((abs(tmp(:,:,1:3:end))));
    colormap parula
    
    % add the squares of the coil images
    im = im + tmp.^2;
    drawnow
    pause(0.5)
end

figure
im = abs(sqrt(im));
%ov([],im, round(dim(1)/2), round(dim(2)/2),round(dim(3)/2),0); 
lightbox(im);
return
