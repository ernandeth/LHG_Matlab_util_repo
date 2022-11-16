function [out3d dcf kernel] = grid3d_lhg(x,y,z,signal, xc, yc, zc, kchoice,  dcf)
% function out3d = grid3d_sinc(x,y,z,signal, xc, yc, zc [, dcf])
%
% 3d gridding interpolation 
%
% x y z    :    are the coordinates of the data that we have
%               dimensions are Nsamples x 1
%
% signal :      is the data that we want to interpolate
%               dimenesions are Nsamples x 1
%
% xc yc zc :    are cartesian coordinates of all the points in the grid
%               dimensions are Nx X Ny X Nz
%
% kchoice:      choice of interpolation kernel 
%               1 = sinc
%               2 = bessel
%               3 = gaussian
%
% dcf :         precalculated density compensation function 
%               use empty matrix if you want to calculate it now
%
% kernel:       interpolation kernel used, but calculated for the dimensions
%               of the kspace data.  Use its FFT to compensate the image
%               intensity in the image domain

outdims = size(xc);

% size of a k-space voxel
dx = abs(yc(2,1,1)-yc(1,1,1));
R = sqrt(xc.^2 + yc.^2 + zc.^2);

%signal = circshift(signal, 2);

% use only the samples that are within V pixels for interpolation
V = 4;
Rmax = V*dx;

dim3 = size(xc);
xc = xc(:);
yc = yc(:);
zc = zc(:);
R = R(:);
out1d = zeros(size(xc));
outcnt = zeros(size(xc));
kernel = zeros(size(xc));

NF = sum(sinc(-V:0.01:V));

outdims = size(xc);
Nsamples = length(signal);
Npix = length(xc);

if isempty(dcf)
    fprintf('\n\tCreating sampling density function ...')
    % inverse of average distance to other samples
    cnt = zeros(Nsamples, 1);
    df = zeros(Nsamples, 1);
    
    parfor n = 1:Nsamples
        
        for nn = 1:Nsamples
            dist = sqrt(...
                (x(n) - x(nn))^2 + ...
                (y(n) - y(nn))^2 + ...
                (z(n) - z(nn))^2);
            
            if dist < Rmax
                df(n) = df(n) + dist ;
                cnt(n) = cnt(n)+1;
            end
        end
    end
    % average distance to neighbors
    dcf = df ./ cnt;
    dcf = dcf-min(dcf);    
    %dcf = dcf.^1.5;
    dcf = 0.999*dcf/max(abs(dcf)) + 0.001;  % normalize to 1, but add a small number to avoid dividing by zero
    
    % Include radius weighting?
    % dcf = (0.01+R/max(R)).^0.5 .* dcf.^3;
    
    fprintf('...DCF done')
else
    fprintf('\n\t(Using pre-computed density function)')
    
end
%
fprintf('\n\tBegin Interpolating ...')


parfor p = 1:Npix

    % convolution kernel - save this for compensation later
    switch kchoice
        case 1  % sinc interpolation kernel
            %H = 0.5*(1-cos(2*pi*dist/dx/2)); % hanning window on the sinc
            kernel(p) = sinc(8*R(p)/dx);
        case 2    % Bessel Function kernel
            kernel(p) = 1/besseli(0, 2.5*R(p)/dx );
        case 3    % Gaussian Function kernel
            kernel(p) = exp(-(R(p)/dx)^2 /4 );
    end
    
    for n=1:Nsamples
        dist = sqrt(...
            (xc(p) - x(n))^2 + ...
            (yc(p) - y(n))^2 + ...
            (zc(p) - z(n))^2);
        
        
        % use only the samples that are within Rmax interpolation
        if dist < Rmax
            
            switch kchoice
                case 1  % sinc interpolation kernel
                    %H = 0.5*(1-cos(2*pi*dist/dx/2)); % hanning window on the sinc
                    W = sinc(8*dist/dx);
                    
                case 2    % Bessel Function kernel 
                    W = 1/besseli(0, 2.5*dist/dx );
                    
                case 3    % Gaussian Function kernel 
                    W = exp(-(dist/dx)^2 /4 );

            end
            out1d(p) = out1d(p) + W*signal(n) .* dcf(n);
            outcnt(p) = outcnt(p) + 1;
            
        end
    end
    %{
    if mod(n,500)==0
        lightbox(log(abs(reshape(out1d./outcnt , dim3)))); 
        title(['n = ' num2str(n) ]);
        drawnow
    end
    %}
end

% average the weighted sums over number of contributing neighbors
out1d = out1d./(outcnt);

outcnt = reshape(outcnt, dim3);
out3d = reshape(out1d, dim3);
kernel = reshape(kernel, dim3);

out3d(isnan(out3d))= eps;
out3d(isinf(out3d))= eps;
kernel(isnan(kernel)) = eps;

%lightbox(log(abs(out3d ))); colormap jet; drawnow
%title('Regridded to Cartesian')

fprintf('\n\t .. done Interpolating. ')

return