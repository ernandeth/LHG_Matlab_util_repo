% Function to Reconstruct Raw data from a multi slice Gems Sequence.
% This function can only reconstruct images sampled on a Cartesian grid.
% [imageMatrix, kSpaceMatrix, procpar] = reconstructGREImage(fileLocation)
%      Inputs:  fileLocation: String. Location of the .fid file containing
%               the FID file of the acquired images
%      Outputs: imageMatrix: NxN complex matrix which is the reconstructed image
%               kSpaceMatrix: NxN complex matrix which is the cartesian samples of the object in K-Space. 
%               procpar:  Structure containing Vnmrj imageing parameters.


function [imageMatrix, kSpaceMatrix, procpar] = reconstructGREImage(fileLocation)
% [imageMatrix, kSpaceMatrix, procpar] = reconstructImage([fileLocation/filename])
% Function to Reconstruct Raw data from a multi slice Gems Sequence
% Cannot yet account for averaging.

% Read in Raw Data
[dta, procpar] = getFid(fileLocation);
 
X = procpar.nv;                  % Number of phase encodes
Y = procpar.np/2;        % Number of points in the readout direction
ns = procpar.ns;  % Number of Slices
kSpaceMatrix = zeros(Y,X,ns);  % Preallocate image matrix
averages = procpar.nt;        % # of averages

% Reshape data into slices (which are along the second dimension)
newData = reshape(dta,[Y,ns,X]);

% Rearrage data so slices are in third dimension and transform into an
% image
for l = 1:ns;
    kSpaceMatrix(:,:,l) = squeeze(newData(:,l,:));
end

imageMatrix = fftshift(ifft2(kSpaceMatrix));
imageMatrix = ifftshift(imageMatrix,3);
end
