function img = ksample2(kxy, object, fov)

% function img = ksample2(kxy, object, FOV)
%
% this function samples an object (2D matrix), given a k-space trajectory
% by computing the signal equation at every voxel.
%
% the trajectory must be given as a vector of complex numbers
% the result is a gridded k-space data set
%
% every pixel in the object is assumed to be 1x1 cm.
%
% the output of this version is a reconstructed image.

GAMMA = 26752;  % rad/s/Gauss

    % the size of the k space data is hard wired and thus the
    % trajectory ius scaled to fit in it.
    matrix_size = size(object,1);
    kmax = 1/fov; % the maximum kspace location is the highest resolution (samples/ cm)
    rstep= fov/size(object,1);  % scaling factor for voxel size
    
    % add the signals from all the voxels
    signal = zeros(size(kxy));
    for xi=1:size(object,1)
        for yi=1:size(object,2)
            
            % compute the contribution to the signal:
            signal = signal + ...
                object(xi, yi) * exp(-i * 2*pi .* ...
                (real(kxy)*(xi-matrix_size/2)*rstep + ...
                 imag(kxy)*(yi-matrix_size/2)*rstep ));
            
        end
    end
    
    
%   an attempt to spped things up.  Actually a lot slower!!!!
%
%     kxy = kxy';
%     % make the object a row vector of voxels
%     object = reshape(object,1, matrix_size*matrix_size);
%     [x,y] = meshgrid([-matrix_size/2: matrix_size/2-1]* fov) ;
%     x = reshape(x,1, matrix_size*matrix_size);
%     y = reshape(y,1, matrix_size*matrix_size);
%     
%     % make the signal a matrix for spiral trajectories ...
%     signal = zeros(max(size(object)), max(size(kxy)));
%         
%     coeffs = exp(-i * 2*pi *( real(kxy)*(x-matrix_size/2) + imag(kxy)*(y-matrix_size/2)));
%     signal = object * coeffs';
%     



    % grid the spiral data and reconstruct
    img = k2image(64*real(kxy)', 64*imag(kxy)', signal.', ones(size(kxy')), 64,2);
    
return