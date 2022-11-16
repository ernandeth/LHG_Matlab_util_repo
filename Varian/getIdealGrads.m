% function  [Gxi, Kxi, Gyi, Kyi, ts, Gx, Kx, Gy, Ky, tg] = getIdealGrads(procpar)
%
%   Function for getting the calculated ideal gradients and their
%   respective K-space trajectories.  The function then goes on to upsample
%   both to the same time step as the sample acquisition.  Gradient
%   waveforms are assumed to be in G/cm
%
% %         Inputs:  procpar:  Procpar file from the spiral scan
%
%           Outputs: G:  Gradient waveforms upsampled by interpolation to
%                        the same grid as (1:sprialLength)/sw
%                    K:  Time integral of Gradient wavform.  K =
%                        gamma*cumsum(G)*gradTimeStep;  interpolated to the same
%                        grid as G.

function  [Gxi, Kxi, Gyi, Kyi, ts, Gx, Kx, Gy, Ky, tg] = getIdealGrads(procpar)
% form is assumed to be in G/cm.

% gamma
gamma = 4257.6;  % Gamma for 1H.
gradTimeStep = procpar.spResolution; % s

%  Get data
[Gx Gy] = getGloverSpirals(procpar);
tg = gradTimeStep*(1:length(Gx));  % Gradient time series

% Account for Spiral in
if strcmp(procpar.spiral_in,'y')
    Gx = -Gx(end:-1:1);
    Gy = -Gy(end:-1:1);
end

% Calculate Kx trajectory
Kx = zeros(size(Gx)); 
Kx(1) = Gx(1)*gamma*gradTimeStep;
for q = 2:length(Gx)
    Kx(q) = Kx(q-1)+Gx(q)*gamma*gradTimeStep;
end

% Calculate Ky trajectory
Ky = zeros(size(Gy)); 
Ky(1) = Gy(1)*gamma*gradTimeStep;
for q = 2:length(Gy)
    Ky(q) = Ky(q-1)+Gy(q)*gamma*gradTimeStep;
end

% Account for Spiral in
sw = procpar.sw;
spiralLength = procpar.np/2;

ts = (1:spiralLength)./sw;

% Correct by summing over the first moment
if strcmp(procpar.spiral_in,'y')
    Kx = Kx-mean(Kx);
    Ky = Ky-mean(Ky);
    ts = ts+abs(max(tg)-max(ts));
end

% Interpolated waveforms
Kxi = interp1(tg,Kx,ts,'spline'); Kxi = Kxi';
Gxi = interp1(tg,Gx,ts,'spline');Gxi = Gxi';
Kyi = interp1(tg,Ky,ts,'spline'); Kyi = Kyi';
Gyi = interp1(tg,Gy,ts,'spline');Gyi = Gyi';


end
