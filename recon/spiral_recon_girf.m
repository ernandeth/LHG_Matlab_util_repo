function  [ pictureCorrected, pictureRaw] = spiral_recon_girf(spiralData, girf_file)
% function  [ pictureCorrected, pictureRaw] = spiral_recon_girf(spiralData, girf_file)

% Set Constants
gamma = 4257.6;  %radians/s/G
dt = 4e-6;     % Gradient sample time in seconds

% Get spiral data.
[dtaArray, procpar] = getSpiralFid( spiralData);

sw = procpar.sw;
ts = (1:size(dtaArray,1))/sw;
dt = ts(2) - ts(1);


%get the ideal gradient waveform--correct with Girf estimate.
%Variables with  i at the end are sampled at 1/sw.  Other variables are
%sampled at dt = 4e-6 seconds.
%ts is time sampled at the same rate as the acquisition bandwidth
% tg is time sampled at the same rate as the gradient amplifiers
% resolution.


[Gyi, Kyi, Gxi, Kxi, ts, Gy, Ky, Gx, Kx, tg] = getIdealGrads(procpar);

KxiHold = Kxi;
KyiHold = Kyi;
GxiHold = Gxi;
GyiHold = Gyi;

q=1;

for delay =0
    
    % Modify Gx and Gy with the GIRF
    load(girf_file)
    
    
    
    % Raw data
    %delay = 5;
    if delay<0
        dta = dtaArray(1:end+delay,:,q);
        Kxi = KxiHold(1-delay:end ,:);
        Kyi = KyiHold(1-delay:end ,:);
        Gxi = GxiHold(1-delay:end ,:);
        Gyi = GyiHold(1-delay:end ,:);
    else
        dta = dtaArray(delay+1:end,:,q);
        Kxi = KxiHold(1:end-delay,:);
        Kyi = KyiHold(1:end-delay,:);
        Gxi = GxiHold(1:end-delay,:);
        Gyi = GyiHold(1:end-delay,:);
    end
    
    % Apply the GIRF
    Gxi_2 = apply_girf(Gxi, xre, xim, xmxf, ts(2)-ts(1));
    Gyi_2 = apply_girf(Gyi, yre, yim, zmxf, ts(2)-ts(1));
    
    Kxi_2 = cumsum(Gxi_2)*(ts(2)-ts(1))*gamma;
    Kyi_2 = cumsum(Gyi_2)*(ts(2)-ts(1))*gamma;
    
    Ki = [];
    Ki_2 = [];
    
    rotang = 2*pi/procpar.Nint;
    rotmat = [ ...
        cos(rotang)   -sin(rotang) ;
        sin(rotang)   cos(rotang)];
    
    tmp1 = [Kxi Kyi];
    tmp2 = [Kxi_2 Kyi_2];
    
    for r=1:procpar.Nint
        
        tmp1 = tmp1 * rotmat ;
        Ki = [Ki; tmp1];
        
        tmp2 = tmp2 * rotmat ;
        Ki_2 = [Ki_2; tmp2];
        
    end
    
    % Regrid Data.  Regrid using both the measured trajectories, and the ideal
    % trajectories corrected by the GIRF.
    
    % Pre-Allocate matrices
    pictureCorrected = zeros(2*procpar.Mat, 2*procpar.Mat, procpar.ns);
    RegridCorrected = pictureCorrected;
    RegridRaw = pictureCorrected;
    pictureRaw = pictureCorrected;
    
    
    
    % delay = -5
    % GIRF corrected trajectory
    %delay =0;
    delay
    
    
    RegridCorrected(:,:,q) = matlabRegrid(Ki_2, (dta), procpar.lro, procpar.Mat, procpar.Mat*2, procpar);
    pictureCorrected(:,:,q) = ifftshift(ifft2((RegridCorrected(:,:,q))));  % Compute the ifft
    
    
    RegridRaw(:,:,q) = matlabRegrid(Ki, (dta), procpar.lro, procpar.Mat,procpar.Mat*2, procpar);
    pictureRaw(:,:,q) = ifftshift(ifft2((RegridRaw(:,:,q))));  % Compute the ifft
    
    
    figure(44)
    
    subplot(121)
    imagesc(abs(pictureCorrected(:,:,q))); colormap hot
    title('Girf corrected')
    axis image
    
    subplot(122)
    imagesc(abs(pictureRaw(:,:,q))); colormap hot
    title('uncorrected')
    axis image
    display('Recon Done');
    
    pause(0.2)
    %end
end




