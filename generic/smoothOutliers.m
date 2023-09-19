function out = smoothOutliers(raw, threshold);
% function out = smoothOutliers(raw, threshold);
%
% Given a 2D data set, smooth out the outliers along the row direction
% This means replacing outliers with the average of their neighbors
% along the row dimension
%
% outliers are defined as being threshold * stddev(row) away from the mean


% allocate space for output:
out = zeros(size(raw));
mu = mean(raw(2:end-1, :), 1);
sigma = std(raw(2:end-1, :), 1);

% go through all the columns
for col=1:size(out,2)
    % identify spikes at each k-space location:
    % first we detrend the data, then we identify the spikes
    % datacolumn = detrend(raw(:,col));
    datacolumn = (raw(:,col));
    outlier_inds = find( abs(datacolumn-mu(col)) > sigma(col)*threshold );

    % now replace the spikes by the mean of their neighbors
    for r=1:length(outlier_inds)
        % Now replace spikes by a mean of its neighbors (in time) -
        % (linear interpolation)
        % make sure we're not overwriting the field map (the first row)!
        if (outlier_inds(r)>1 &&   outlier_inds(r) < length(datacolumn) -1)

            row = outlier_inds(r);
            out(row, col) = 0.5*(raw(row-1, col) + raw(row+1, col)) ;
        
        end
    end
end

return

