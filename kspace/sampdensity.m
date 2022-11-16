function dens = sampdensity(samplocs,N)
% function dens = sampdensity(samplocs, N)
%
% compute the sampling density of a set of points (samplocs)
% by summing the distance of the points N nearest neighbors
%
% samplocs are the coordinates of the samples usually N x 3
% N is the number of neighbors to evaluate
%
% the output is scaled so that the maximum density is 1
%

fprintf('\nMaking Density compensation function...');

% first downsample the number of locations
F = 10;
N0 = [1:length(samplocs)]';
N1 = [1:F:length(samplocs)]';
samplocs = samplocs(1:F:end,:);
N = floor(N/F);

% Initialize output:
dens = zeros(size(samplocs,1), 1);

for n = 1:length(dens)
    distances = vecnorm((samplocs(n,:) - samplocs), 2, 2);
    distances = sort(distances);
    dens(n) = sum(distances(1:N).^2);
end
dens = dens - dens(min(~isinf(dens)));
dens = dens./max(dens(~isinf(dens)));
dens(isnan(dens)) = 0;
dens(isinf(dens)) = 1;

% upsample back to the original number
% we assume that the desity function is smooth
dens = abs(interp1(N1,dens,N0));
return