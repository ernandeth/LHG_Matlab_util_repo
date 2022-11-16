%    Function using matlab's griddata command to regrid spiral data onto a
%    cartesian grid.  
%   
% function  kSpace = matlabRegrid(K, spData, FOV, Mat, gridSize, procpar)
%           Inputs:  K: complex Nx1 or Nx2 vector in the of each sample's K-space location
%                       IE:  K = [Kx Ky] or Kx + i*Ky.  Units: cycles/cm.
%                    spData:  Nx1 complex points of the estimate of K space at the
%                              points in K.  
%                    FOV:  Field of view of the acquistion Units: cm
%                    Mat:  Intended size of the reconstructed image.
%                    gridSize: size of the grid on which to interpolate the reconstructed data.
%                    procpar:  Vnmrj paramter file
% 
%           Outputs:  kSpace: Interpolation of spiralData onto a gridSize X gridSize array of complex points
%                             Gridsize covers locations 
%                             linspace(-Mat/FOV/2,Mat/FOV/2,gridSize);


function  kSpace = matlabRegrid(varargin)
% function  kSpace = matlabRegrid(K,spData,FOV,Mat,gridSize, procpar)

% Formatting and parsing
if  numel(varargin) < 5
    display('kSpace = matlabRegrid(kSpaceLocations,spiralData,FOV,Mat,gridSize)');
    error('Not enough arguments -- 5 arguments required. ');
end

K  = varargin{1};
spData = varargin{2};
FOV = varargin{3};
Mat = varargin{4};
gridSize = varargin{5};
procpar = varargin{6};

% format K into Nx2 vector
[m n] = size(K);
if n == 1
    U = K;
    K(:,1) = real(U);
    K(:,2) = imag(U);
elseif n >= 3
    error('K space data must either be complex or Nx2. IE: [xcoords, ycoords]') 
else
end

[m,n] = size(spData);
if n >=2
    spData = reshape(spData,numel(spData),1);
end

% Form a grid on which to grid the spiral data
KXi = linspace(-Mat/FOV/2,Mat/FOV/2,gridSize);
[KXi,KYi] = meshgrid(KXi,KXi);

% Use griddata to interpolate the scattered data
kSpace1 = griddata(K(:,1)',K(:,2),spData,KXi,KYi,'cubic'); 
kSpace1(isnan(kSpace1)) = 0;

% Add a linear shift due to the imaging window offset
y0 = procpar.pro;
x0 = procpar.ppe;
fx = 1i*2*pi*x0*KXi;
fy = 1i*2*pi*y0*KYi;
kSpace = kSpace1; %.*exp(fx).*exp(fy);

end