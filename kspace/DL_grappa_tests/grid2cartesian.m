function cart_signal = grid2cartesian(kx, ky, kz, data, dim3, ROI_fraction)
% function cart_signal = grid2cartesian(kx, ky, kz, data, dim3, ROI_fraction)
%
% function to interpolate multi-channel 3D data into a cartesian sampling grid
% from an arbitrary sampling pattern specified by kx, ky, kz
% inputs:
%
%   kx, ky, kz: locations of the data in the sampled grid
%   dim3:  is the sizes of the output data grid
%   data:  Nsamples x Nchannels
%   ROI_fraction : this is a fraction of the data that you want to
%       interpolate.  If you want the whole grid, use 1, if you just want the
%       center half, use 0.5
%
% output:
%   cart_signal : 3D cartesian sampled grid, arranged in columns
%       dim3*dim3*dim3  x Nchannels
%

% figure out how many coisl are present
Ncoils = size(data,2);


fprintf('\rGRIDDING: Making cartesian calibration data from non cartesian data...\n')
% now interpolate spiral data into a cartesian grid.  The GRAPPA
% interpolation will be done on these data.
% we have to do this for each coil

% create the mesh for interpolation
cart_signal = [];
dcf = [];
ROI_dim = floor(dim3*ROI_fraction);    % size of cartesian calibration data matrix   
R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
ROI_fraction = ROI_fraction * Kmax;    % convert from fraction to k space units.

[subkx , subky, subkz ]= meshgrid( ...
    linspace(-ROI_fraction,ROI_fraction,ROI_dim(1)), ...
    linspace(-ROI_fraction,ROI_fraction,ROI_dim(2)), ...
    linspace(-ROI_fraction,ROI_fraction,ROI_dim(3)));
im = zeros(size(subkx));

% Make a Gaussian filter in k-space.  WIll use this for 
% smoothness as regularization.
% Gfilter = exp(-(subkx.^2 + subky.^2 + subkz.^2) ./ ROI_dim(1).^2 );

for coilnum=1:Ncoils
     
    % calibration region for each coil
    fprintf('\rcoil  %d ... ', coilnum);
    
    %{
    tic
    % grid the calibration dat into a cartesian grid
    [tmp dcf kernel] = grid3d_lhg(...
        kx, ky, kz, ...
        data(:,coilnum), ...
        subkx, subky, subkz, ...
        2, dcf);
    toc
    %}
    
    % Alternatively ... use griddata
    tmp = griddata(kx, ky, kz, data(:,coilnum), subkx, subky, subkz,'linear');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    
    % apply a filter to the data -  smoothness in image domain
    % this should act as a regularizer of sorts.
    % tmp = tmp .* Gfilter;
    % tmp = smooth3(tmp, 'gaussian',5);  <---- this never seems to help

    % testing images for recon
    im = im + (fftshift(abs(ifftn(tmp)))).^2;
    %fk = fftshift(abs(fftn(kernel)));
    
    subplot(2,1,1)
    lightbox(log(abs(tmp ))); 
    title(['coil ' num2str(coilnum) ' calibration region (gridded)']); 
    
    subplot(2,1,2)
    lightbox(sqrt(im)); 
    title(['coil: ' num2str(coilnum) ' cal region in image domain']); 
    
%     subplot(3,1,3)
%     %lightbox(abs(fk(10:end-10,10:end-10, 10:end-10)),[],5,[]); 
%     title('kernel rolloff (image domain)');
    
    drawnow
    
    cart_signal = [cart_signal tmp(:)];
end

return