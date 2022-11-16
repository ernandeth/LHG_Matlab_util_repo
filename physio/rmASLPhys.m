function PM = rmASLPhys(rootname, threshold, regs, maskFile);
% function rmASLPhys(rootname , threshold ,regs, maskFile);
% 
% identify CSF voxels as the brightest in the image (mean + 2 stadndard
% deviations)
% extract the mean time course in those voxels
% append a baseline to it and a an oscillating baseline to it
%
% threshold :  number of std. devs used to determine
%              which pixels are CSF
% output:      residual.img

if threshold==[]
	threshold = 2.5;
end

fprintf('\nReading data ...');
[data ,h] = read_img_series(rootname);
output = zeros(size(data));

% split the data into control and tag
cdata = data(1:2:end,:);
tdata = data(2:2:end,:);

if isempty(maskFile) 
	% compute the mask based on the control images
	fprintf('\nSlecting voxels at %f SD threshold', threshold);
	mdata = mean(cdata,1);
	md = mean(mdata);
	sd = std(mdata);
	threshold = md + threshold*sd;
	inds = find(mdata(1:h.xdim*h.ydim*h.zdim/2)> threshold);
	mask = zeros(size(mdata));
	mask(inds) = 1;
	%write the results to file for diagnostic purposes
	hm = h;
	hm.tdim = 1;
	fprintf('\nwriting mask, and mean  files ...');
	write_img('mask.img',mask,hm);
	write_img('mean.img',md,hm);
else
	mask = read_img(maskFile);
	inds = find(mask);
end


h.tdim = size(data,1);
if isempty(regs)
	% extract the time series from the mask
    % for making the Physio regressors
	ts = data(:,inds);
	ts = mean(ts,2);
	% build design matrix for regression
	% baseline, baseline ASL, linear trend, CSF timecourse
	PM = ones(h.tdim,4);
	PM(2:2:end,2) = -1;
    PM(:,3) = linspace(-1, 1, length(ts));
	PM(:,4) = ts - mean(ts);
    	% derivative regressor
	% dts = diff(ts); dts = dts - mean(dts);
	% PM(2:end,4) = dts;  PM(1,4) = 0;
    	% linear trend regressor
    	% PM(:,4) = linspace(-1, 1, length(ts))
	contrast = [0 0 1 1];
else
	PM = ones(h.tdim,4);
	%PM(2:2:end,2) = -1;
	%PM(:,3:end) = regs;
	PM(:,2:end) = regs;
	contrast = ones(1,4);
end
contrast = repmat(contrast,h.tdim,1);
% estimate models
beta_est = pinv(PM) * data;
% subtract fitted model from data to see how well it fits
fit = PM*beta_est;
residuals= data - fit;
% subtract fitted model from data , except for the first two regressor
beta_est(1:2,:)=0;
output = data - PM(:, 4) * beta_est(4,:);

if 0
for c=1:size(data,2)
	plot(output(:,c),'r')
	hold on
	plot(data(:,c))
	title(sprintf('before: %0.2f  after: %0.2f', std(data(:,c)), std(output(:,c)))  )
	hold off
	drawnow
	pause
end
end
warning off
h.tdim = size(output,1);
write_img(['clean_' rootname '.img'], output, h);
h.tdim = 1;
write_img('residualSTD.img',100*std(residuals,[],1),h)
write_img('originalSTD.img',100*std(data,[],1),h)
write_img('correctedSTD.img',100*std(output,[],1),h)
save rmASLPhysvars beta_est PM
warning on
return
