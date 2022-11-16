% Function to Reconstruct Raw data from a multi slice Gems Sequence acquired with an arrayed parameter.
% For example, acquisitions were made at several echo times or band widths.
% This function can only reconstruct images sampled on a Cartesian grid.
% [imageMatrix, kSpaceMatrix, procpar] = reconstructGREImageArrayed(fileLocation)
%      Inputs:  fileLocation: String. Location of the .fid file containing
%               the FID file of the acquired images
%      Outputs: imageMatrix: NxN complex matrix which is the reconstructed image
%               kSpaceMatrix: NxN complex matrix which is the cartesian samples of the object in K-Space. 
%               procpar:  Structure containing Vnmrj imageing parameters.



function [imageMatrix, kSpaceMatrix, procpar] = reconstructGREImageArrayed(fileLocation)
% [imageMatrix, kSpaceMatrix, procpar] = reconstructImage(fileLocation)
% Function to Reconstruct Raw data from a multi slice Gems Sequence
% Cannot yet account for averaging.

% Read in Raw Data
cd(fileLocation);  % Move to location of fid file
[dta, fheader, bheader]=readfid('fid',-1);  %Read in fid file
procpar = tryreadprocpar();   % Read in image parameters
 
X = procpar.nv;                  % Number of phase encodes
Y = fheader.np/2;        % Number of points in the readout direction
ns = procpar.ns;  % Number of Slices
kSpaceMatrix = zeros(Y,X,ns);  % Preallocate image matrix
averages = procpar.nt;        % # of averages

% Reshape data into slices (which are along the second dimension)
% newData = reshape(dta,[Y,ns,X]);
kSpaceMatrix = reshape(dta,[Y,X,2]);

% Rearrage data so slices are in third dimension and transform into an
% image
% for l = 1:ns;
%     kSpaceMatrix(:,:,l) = squeeze(newData(:,l,:));
% end

imageMatrix = fftshift(ifft2(kSpaceMatrix));
imageMatrix = ifftshift(imageMatrix,3);
end
