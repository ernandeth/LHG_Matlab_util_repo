function dcf = DCF_vsasl3dflex(ks,dcfP,dcfN,doRecycle,numWorkers)
% function dcf = DCF_vsasl3dflex(ks,dcfP,dcfN,doRecycle,numWorkers)
    
    % Get dimensions of trajectory
    ndat = size(ks,1);
    nleaves = size(ks,3);
    nslices = size(ks,4);
    
    % Initialize 'last' parameters for if statements
    lastdcf = []; lastdcfN = [];
    
    % Define exponential distance weighting function
    rfun = @(x) x.^dcfP;
    
    % If doRecycle flag is on & the file exists, load the previous dcf
    %   calculation to save computation
    if doRecycle && isfile('lastdcf.mat')
        load lastdcf.mat lastdcf lastdcfN
    end
    
    % Check that dimensions and parameters of previous dcf are consistent
    if isequal(size(lastdcf),[ndat,nleaves,nslices]) && isequal(dcfN,lastdcfN)
    
        % ... if so, set dcf equal to previous dcf
        dcf = lastdcf;
        fprintf('\nSuccessfully loaded previous density compensation function');
        fprintf(' with dcfN = %d',dcfN);
    
    else % Otherwise recalculate it

        % Initialize dcf output
        dcf = zeros(ndat,nleaves,nslices);

        % Set up progress message
        fprintf('\nCreating sampling density compensation function');
        fprintf(' with dcfP = %.2f and dcfN = %d, ',dcfP,dcfN);
        str = sprintf('leaf %d/%d, slice %d/%d...',1,nleaves,1,nslices);
        fprintf(str);

        % Index each sample
        for leafn = 1:nleaves
            for slicen = 1:nslices

                % Print progress message
                spiraln = (leafn-1)*nslices + slicen;
                if spiraln>1
                    fprintf(repmat('\b',1,length(str)));
                    str = sprintf('leaf %d/%d, slice %d/%d...',leafn,nleaves,slicen,nslices);
                    fprintf(str);
                end

                % Calculate density at each sample
                parfor (datn = 1:ndat, numWorkers)
                    % Calculate euclidean distance between sample and all
                    % other samples
                    dists = vecnorm(ks(:,:,:,:) - ks(datn,:,leafn,slicen), 2, 2);

                    % Set dcf at sample as average distance to closest dcfN
                    % neighbors
                    dcf(datn,leafn,slicen) = mean(mink(dists(:),dcfN));
                end

            end
        end
        
        dcf(isinf(dcf)) = max(dcf(~isinf(dcf)),[],'all');
        dcf(isnan(dcf)) = min(dcf(~isnan(dcf)),[],'all');
    
    end
    
    % Save to file for recycling
    if doRecycle
        lastdcf = dcf;
        lastdcfN = dcfN;
        save lastdcf.mat lastdcf lastdcfN
    end
    
    % Normalize dcf, apply radial function, and remove inf/nan values
    dcf = rfun((dcf - min(dcf,[],'all'))/(max(dcf,[],'all') - min(dcf,[],'all')) * 0.95 + 0.05);
    dcf(isinf(dcf)) = 1; dcf(isnan(dcf)) = 0;    

end