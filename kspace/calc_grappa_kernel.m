function [GWeights , Rnk, RMSE] = calc_grappa_kernel(TargetPoints, Neighbors)
% function [GWeights , rank, RMSE]  = cal_ncgrappa_kernel(TargetPoints, Neighbors)
%
% solve inverse problem for interpolation kernel:
%   TargetnPoints = Neighbors*GWeights
% so ...
%   GWeights = inv(Neighbors) * TargetPoints ;
%  we use a Pseudo INV using regularization lambda:
%
%   neighbors    % dimensions: patchsize x Ncoils*Nnbrs
%   Weights      % dimensions: Ncoils*Nnbrs x Ncoils
%   targets      % dimensions: patchsize x Ncoils
%
% return the weights: Ncoils*Nnbrs


Nnbrs = size(Neighbors,1);
Ncoils = size(Neighbors,2);

lambda = 1e-6;
tol = 1e-18;
A = Neighbors;
B = TargetPoints;

%GWeights = (inv(A'*A + lambda*eye(size(A,2)))*A') * B;
% GWeights = pinv(A, tol)*B;
GWeights = pinv(A)*B;

% diagnostics:
Rnk= rank(Neighbors);

RMSE = (B - A*GWeights) ./ B ;
RMSE = abs(RMSE);
RMSE = mean(RMSE.^2, 1 , 'omitnan');
RMSE = 100*sqrt(RMSE);
return