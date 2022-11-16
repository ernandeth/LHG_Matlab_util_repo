function  [pictureMeasured, pictureCorrected, pictureRaw] = reconArraySpiralsWith_Without_GIRF(varargin)
% function  [pictureMeasured, pictureCorrected, pictureRaw] = reconArraySpiralsWith_Without_GIRF(varargin)



% Parse inputs
td = varargin{1};               % Top directory
spiralData = varargin{2};
if numel(varargin)>=3
    KzDir = varargin{3};
    KxDir = varargin{4};
    girf_file = varargin{5};
end

% Set Constants
gamma = 4257.6;  %radians/s/G
dt = 4e-6;     % Gradient sample time in seconds

% Get spiral data.
[dtaArray, procpar] = getSpiralFidArray([td '/' spiralData]);

sw = procpar.sw;
ts = (1:size(dtaArray,1))/sw;
dt = ts(2) - ts(1);

% get Kspace locations and upsample to correct number of points if
% necessary
if numel(varargin) >=3;
    KzMeas = getTraj([td '/' KzDir]);
    KxMeas = getTraj([td '/' KxDir]);
    
    KzMeas = -KzMeas;
    KxMeas = -KxMeas;
end

%get the ideal gradient waveform--correct with Girf estimate.
%Variables with  i at the end are sampled at 1/sw.  Other variables are
%sampled at dt = 4e-6 seconds.
%ts is time sampled at the same rate as the acquisition bandwidth
% tg is time sampled at the same rate as the gradient amplifiers
% resolution.

%[Gxi, Kxi, Gzi, Kzi, ts, Gx, Kx, Gz, Kz, tg] = getIdealGrads(procpar);
% Note the swapping of gradient axes!

[Gzi, Kzi, Gxi, Kxi, ts, Gz, Kz, Gx, Kx, tg] = getIdealGrads(procpar);

% Modify Gx and Gz with the GIRF
load(girf_file)

Gxi_2 = apply_girf(Gxi, xre, xim, xmxf, ts(2)-ts(1));
Gzi_2 = apply_girf(Gzi, zre, zim, zmxf, ts(2)-ts(1));

Kxi_2 = cumsum(Gxi_2)*(ts(2)-ts(1))*gamma;
Kzi_2 = cumsum(Gzi_2)*(ts(2)-ts(1))*gamma;


% Rotate trajectories for multishot imaging
[Ki] = rotTraj([Kxi_2,Kzi_2],procpar.Nint);
[Gi] = rotTraj([Gxi_2,Gzi_2],procpar.Nint);

if numel(varargin)>=3
    [KMeas] = rotTraj([KxMeas, KzMeas],procpar.Nint);
end

% Compute the transfer function comparing the nominal and measured gradient
% waveforms.
GxMeas = (KxMeas(2:end) - KxMeas(1:end-1))/dt/gamma;
Hx = fft(GxMeas) ./ fft(Gxi(2:end));
GzMeas = (KzMeas(2:end) - KzMeas(1:end-1))/dt/gamma;
Hz = fft(GzMeas) ./ fft(Gzi(2:end));


%{
figure(31)
subplot(211)
plot(Kxi);
hold on
plot(Kxi_2,'r');
plot(KxMeas,'g');
legend('Nominal', 'GIRF Corrected', 'Measured','Location', 'SouthWest');
hold off

subplot(212)
plot(Kzi);
hold on
plot(Kzi_2,'r');
plot(KzMeas,'g');
legend('Nominal', 'GIRF Corrected', 'Measured','Location', 'SouthWest');
hold off
%}

% Regrid Data.  Regrid using both the measured trajectories, and the ideal
% trajectories corrected by the GIRF.

% Pre-Allocate matrices
RegridMeasured = zeros(2*procpar.Mat, 2*procpar.Mat, 1);  % Preallocate Memory
pictureMeasured = zeros(size(RegridMeasured));
RegridCorrected = RegridMeasured;
pictureCorrected = pictureMeasured;
RegridRaw = RegridCorrected;
pictureRaw = pictureMeasured;

%for q = 1:1
KiHold = Ki;
KMeasHold = KMeas;
KxiHold = Kxi;
KziHold = Kzi;

for delay = -55: 10:-5
    q=5;
    % delay = -5
    % GIRF corrected trajectory
    %delay =0;
    delay
    if delay<0
        dta = dtaArray(1:end+delay,q);
        Ki = KiHold(1-delay:end ,:);
    else
        dta = dtaArray(delay+1:end,q);
        Ki = KiHold(1:end-delay ,:);
    end
    RegridCorrected(:,:,q) = matlabRegrid(Ki,(dta),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    pictureCorrected(:,:,q) = ifftshift(ifft2((RegridCorrected(:,:,q))));  % Compute the ifft
    
 
    % Raw data
    delay = 6;
   if delay<0
        dta = dtaArray(1:end+delay,q);
        Kxi = KxiHold(1-delay:end ,:);
        Kzi = KziHold(1-delay:end ,:);
    else
        dta = dtaArray(delay+1:end,q);
        Kxi = KxiHold(1:end-delay,:);
        Kzi = KziHold(1:end-delay,:);
    end
 
    RegridRaw(:,:,q) = matlabRegrid([Kxi, Kzi],(dta),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    pictureRaw(:,:,q) = ifftshift(ifft2((RegridRaw(:,:,q))));  % Compute the ifft
    
        % measured trajectury
    delay = 0;
    dta = dtaArray(:,q);
    RegridMeasured(:,:,q) = matlabRegrid(KMeas,(dta),procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    pictureMeasured(:,:,q) = ifftshift(ifft2((RegridMeasured(:,:,q))));  % Compute the ifft
    
    figure(44)
    subplot(131)
    imagesc(abs(pictureMeasured(:,:,q))); colormap hot
    title('Measured')
    axis image
    
    subplot(132)
    imagesc(abs(pictureCorrected(:,:,q))); colormap hot
    title('Girf corrected')
    axis image
    
    subplot(133)
    imagesc(abs(pictureRaw(:,:,q))); colormap hot
    title('uncorrected')
    axis image
    display('Recon Done');
    
    pause(0.2)
%end
end




