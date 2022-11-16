%   Function for getting the calculated ideal gradients and their
%   respective K-space trajectories.  
%   The function then goes on to modulate the frequency content of the gradients to correct 
%   for the impulse response function of the gradient hardware (estimated elsewhere).
%   The function then goes on to upsample
%   both to the same time step as the sample acquisition.  Gradient
%   waveforms are assumed to be in G/cm
%
% [G, K] = getIdealGrads(fileLocation, fileName, gradTimeStep, sw, spiralLength);  
% %         Inputs:  procpar:  Procpar file from the spiral scan
%                    
%           Outputs: G:  Gradient waveforms upsampled by interpolation to
%                        the same grid as (1:sprialLength)/sw
%                    K:  Time integral of Gradient wavform.  K =
%                        gamma*cumsum(G)*gradTimeStep;  interpolated to the same
%                        grid as G.

function  [Gxi, Kxi, Gyi, Kyi, ts, Gx, Kx, Gy, Ky, tg] = getIdealGradsCorrectWithGIRF(procpar)
% function  [G, K] = getIdealGrads(fileLocation, fileName,gradTimeStep, sw, spiralLength)  Gradient wave
% form is assumed to be in G/cm.

% gamma
gamma = 4257.707747;  % Gamma for 1H.
gradTimeStep = procpar.spResolution; % s

%  Get data
[Gx1 Gy1] = getGloverSpirals(procpar);
tg = gradTimeStep*(1:length(Gx1));  % Gradient time series

% Apply GIRF Correction

% Load in GIRF Parameters
load(['C:\Users\MR-Histo\Documents\My Box Files\Research\Projects'...
       '\Spiral on 7T\Gradient Impulse Response Correction\girfParameters.mat']);
   
%    FFT Gx and Gy
GX = fftshift(fft(Gx1));
GY = fftshift(fft(Gy1));
GY1 = GY; GX1 = GX;
df = 1/(tg(end));
Nyq = (1/gradTimeStep/2);
frq = (-Nyq):df:(Nyq)-df;        % Spectrum of Gx and Gy.
frqConstruct = min(abs(frq)):df:cutoff;      % Spectrum to reconstruct the transfer function

% Find the indices that match the range of the GIRF
idx1  = find(frq<=cutoff & frq>= -cutoff);

% Compute the GIRF over the desired frequency range
pRealX = allrealp(3,:);
pImagX = allimagp(3,:);
pRealY = allrealp(2,:);
pImagY = allimagp(2,:);

H_EstimatedXReal = polyval(pRealX,frqConstruct);
H_EstimatedXImag = polyval(pImagX,frqConstruct);
H_EstimatedYReal = polyval(pRealY,frqConstruct);
H_EstimatedYImag = polyval(pImagY,frqConstruct);

H_EstimatedXReal = [H_EstimatedXReal(end:-1:1),H_EstimatedXReal];
H_EstimatedXImag = [-H_EstimatedXImag(end:-1:1),-H_EstimatedXImag];
H_EstimatedYReal = [H_EstimatedYReal(end:-1:1),H_EstimatedYReal];
H_EstimatedYImag = [-H_EstimatedYImag(end:-1:1),-H_EstimatedYImag];

% Flip over the 

% Apply the GIRF to Gx and Gy
GXReal = real(GX(idx1)).*-H_EstimatedXReal';
GXImag = imag(GX(idx1)).*-H_EstimatedXImag';
GYReal = real(GY(idx1)).*-H_EstimatedYReal';
GYImag = imag(GY(idx1)).*-H_EstimatedYImag';
  
% Convert back to complex
GX(idx1) = GXReal+1i*GXImag;
GY(idx1) = GYReal+1i*GYImag;

% Reduce everything out of this range to zero
GX(~idx1) = 0;
GY(~idx1) = 0;

% Convert back into the time domain

Gx = real(ifft(ifftshift(GX)));
Gy = real(ifft(ifftshift(GY)));



% Account for Spiral in
if strcmp(procpar.spiral_in,'y')
Gx = -Gx(end:-1:1);
Gy = -Gy(end:-1:1);
end

% Calculate Kx trajectory
Kx = zeros(size(Gx)); Kx(1) = Gx(1)*gamma*gradTimeStep;
for q = 2:length(Gx)
    Kx(q) = Kx(q-1)+Gx(q)*gamma*gradTimeStep;
end

% Calculate Ky trajectory
Ky = zeros(size(Gy)); Ky(1) = Gy(1)*gamma*gradTimeStep;
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
