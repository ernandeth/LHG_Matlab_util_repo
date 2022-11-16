function slope = regression(x,y)
% function slope= regression(x,y)
% this function returns the estimate of the betas
% this could be slope and intercept

	x = [ones(size(x,1) , 1)  x];
	c=[0 1]';

	p = size(x,2);
	n = size(x,1);

	xtx_inv = pinv(x'*x);
	beta_est = xtx_inv*x'*y;
	slope=beta_est(2);

	RSS = y - x*beta_est;
	RSS = sum(RSS.^2);

	var_est = RSS/(n-p);
	cov_beta = var_est * xtx_inv;

	% ---- %

	v = var_est*(c' * xtx_inv * c);
	t = (beta_est' * c) ./ sqrt(v);

	save glm.mat beta_est var_est

	% For testing only.  Change this when you are done !!!
	%t = beta_est' * c;
	%t = v;
	%t = RSS;

return
