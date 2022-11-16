function [out3d dcf] = grid3d_sinc(x,y,z,signal, xc, yc, zc, dcf)
% function out3d = grid3d_sinc(x,y,z,signal, xc, yc, zc [, dcf])
%
% 3d gridding interpolation
%
% dcf :         precalculated density compensation function [optional]
%
% xc yc zc :    are cartesian coordinates of all the points in the grid
%               dimensions are Nx X Ny X Nz
%
% x y z    :    are the coordinates of the data that we have
%               dimensions are Nsamples x 1
%
% signal :      is the data that we want to interpolate
%               dimenesions are Nsamples x 1
%

outdims = size(xc);

% size of a k-space voxel
dx = abs(yc(2,1,1)-yc(1,1,1));
R = sqrt(x.^2 + y.^2 + z.^2);

% use only the samples that are within 5 pixels for interpolation
V = 4;
Rmax = V*dx;
dim3 = size(xc);
xc = xc(:);
yc = yc(:);
zc = zc(:);
out1d = zeros(size(xc));
outcnt = zeros(size(xc));

NF = sum(sinc(-V:0.01:V));

outdims = size(xc);
Nsamples = length(signal);
Npix = length(xc);

if nargin < 8
    fprintf('creating sampling density function ...')
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
    % inverse of average distance to other samples
    dcf = cnt;%./df;
    fprintf('...DCF done')
else
    fprintf('(Using pre-computed density function)')
    
end
%
fprintf('\n Begin Interpolating ...')

for n=1:Nsamples
    
    for p = 1:Npix
        dist = sqrt(...
            (xc(p) - x(n))^2 + ...
            (yc(p) - y(n))^2 + ...
            (zc(p) - z(n))^2);
        
        
        % use only the samples that are within Rmax interpolation
        if dist < Rmax
            %             H = 0.5*(1-cos(2*pi*dist/Rmax)); % hanning window on the sinc
            %             W = H*sinc(dist/dx);
            % KB kernel version
            W = 1/besseli(0, dist/dx );
            out1d(p) = out1d(p) + W*signal(n)/dcf(n);
            outcnt(p) = outcnt(p) + 1;
            
        end
    end
    %{
    if mod(n,5000)==0
        lightbox(log(abs(reshape(out1d./outcnt , dim3)))); drawnow
    end
    %}
end
fprintf('\n .. done Interpolating')

% average the weighted sums over number of contributing neighbors
out1d = out1d./(outcnt);

outcnt = reshape(outcnt, dim3);
out3d = reshape(out1d, dim3);
out3d(isnan(out3d))= 0;
lightbox(log(abs(out3d ))); colormap jet; drawnow
title('Regridded to Cartesian')


return