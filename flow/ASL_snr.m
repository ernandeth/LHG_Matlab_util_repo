function [tSNR,sSNR,msk] = ASL_snr(sub,msk,show)

    if nargin<1 || isempty(sub)
        % if subtraction data is not defined, read it from file
        sub = nii2matrix('./sub.nii');
    end

    if nargin<2 || isempty(msk)
        % if msk is not defined, define it as a threshold of 1
        msk = 1;
    end

    if isscalar(msk)
        % if msk is a scalar, create a mask using mean_con and 'msk' as a
        %   threshold
        meancon = nii2matrix('./mean_con.nii');
        msk = 1*imfill(meancon.^2>=msk*std(meancon.^2,[],'all'),'holes');
    end
    
    if nargin<3 || isempty(show)
        % if show is not defined, set it to 1 (default is to show figures)
        show = 1;
    end
    
    % calculate signal and noise values
    tms = mean(sub,4); % temporal mean
    tsds = std(sub,[],4) + eps; % temporal std
    signal = mean(tms(msk==1),'all'); % spatial mean of temporal mean inside the ROI
    tNoise = mean(tsds(msk~=1),'all'); % spatial mean of temporal std outside the ROI
    sNoise = std(tms(msk~=1),[],'all'); % spatial std of temporal mean outside the ROI
    
    % calculate ratios
    tSNR = signal/tNoise;
    sSNR = signal/sNoise;
    
    % show figures
    if show
        
        subplot(2,2,1)
        lightbox(tms.*msk);
        title('Temporal mean inside ROI');
        
        subplot(2,2,2)
        lightbox(tsds.*(~msk));
        title('Temporal std outside ROI');
        
        subplot(2,2,3)
        lightbox(tms.*(~msk));
        title('Temporal mean outside ROI');
        
        subplot(2,2,4)
        lightbox(tms./tsds);
        title('tSNR per voxel (log)');
        
        sgtitle(sprintf('tSNR = %.2f, sSNR = %.2f',tSNR,sSNR));
        
    end
    
end