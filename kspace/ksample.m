function kspace = ksample(kxy, object, fov)

% function kspace = ksample(kxy, object, FOV)
%
% this function samples an object (2D matrix), given a k-space trajectory
% the trajectory must be given as a vector of complex numbers
% the result is a gridded k-space data set
%
GAMMA = 26752;  % rad/s/Gauss

    % the size of the k space data is hard wired and thus the
    % trajectory ius scaled to fit in it.
    matrix_size = size(object,1);
    kmax = 1/fov; %(samples/ cm)
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



    % grid the spiral data
    [kx,ky] = meshgrid([-matrix_size/2: matrix_size/2-1]* (1/fov) );
    %keyboard
    kspace=griddata(real(kxy), imag(kxy), signal, kx, ky); 
    kspace(find(isnan(kspace))) = 0;
    
return