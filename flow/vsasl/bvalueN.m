%% Calculate b-value of a discrete gradient waveform
%
% b = bvalueN(G, GUP, gamrad)
%
% in:
%	   G      - gradient waveform (units/cm)
%	   GUP    - gradient update time (µs)
%      gamrad - gyromagnetic ratio (rad/s/units)
%      output - Output type: 'end' = final b-value, 'all' = b-value
%               time-course with same resolution as G
%
% out:
%      b - calculated diffusion b-value (s/mm^2)
%
% Reference: Bernstein et al., Handbook of MRI Pulse Sequences (2004) page 278
%
% Written by Joseph G. Woods, CFMRI, UCSD, Aug 2020

function b = bvalueN(G, GUP, gamrad, output)

if ~exist('output','var') || isempty(output)
    output = 'end';
end

sizG = size(G);                 % save size of G for output
G    = reshape(G, sizG(1), []); % reshape to have size Nt x []

G_mm  = G / 10;     % convert from units/cm to units/mm
GUP_s = GUP * 1e-6; % convert from µs to s

% Create local sum function based on the output type
switch output
    case 'end'; suml = @(x) sum(x,1);    % sum along time dimension
    case 'all'; suml = @(x) cumsum(x,1); % cumsum along time dimension
end

% Calculate the b-value
b = gamrad^2 * suml( GUP_s * (GUP_s * cumsum(G_mm)).^2 );

% Reshape to have same dimensions as size(G,2:end)
b = reshape(b, [size(b,1), sizG(2:end)] );

end
