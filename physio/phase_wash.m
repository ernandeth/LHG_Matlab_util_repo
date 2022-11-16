function phase_wash(root, rho_thres);
% function phase_wash(root, rho_thres);
%
% this function regresses out the phase from the magnitude in a time series
% of images
%
% it removes the phase of those pixels, whose phase correlates with the
% magnitude above the threshold.
%

[mag h]= read_img(root);
[ph h] = read_img(['p_' root]);
ph = ph/1000;
% num  = diag(mag'*ph) ;
% den = sqrt(sum(mag.^2, 1)) .* sqrt(sum( ph.^2,1));
% rho = num ./ den';
%
% rho = abs(rho);
warning off

Nslices = h.dim(4);
Npix = h.dim(2) * h.dim(3);
Nframes = size(mag,1);
% extract a phase time course for the slice
rho = zeros(1,Npix);
mag2 = zeros(size(mag));
	
for s=1:Nslices
	tmp = zeros(Nframes,1);
	% compute correlation coefficient
	for p = (s-1)*Npix+1: s*Npix

		rho1 = corrcoef(mag(:,p), ph(:,p));
		rho(p) = abs(rho1(1,2));

		if rho(p) > rho_thres
			tmp = tmp + ph(:,p);	
		end
	end

	tmp = tmp - mean(tmp);
	plot(tmp);
	% orthogonalize the labeling waveform from the phase
	reg1 = ones(length(tmp),1);
	reg2 = reg1;
	reg1(1:2:end) = -1;
	reg3 = 1:length(reg1);
	reg3 = reg3'-mean(reg3);

	c = [1 0 0 ; 0 0 0; 0 0 0];
	X = [reg1 reg2 reg3];

	bhat = pinv(X)*tmp;
	tmp = tmp - X *c* bhat;
	hold on ; plot(tmp,'r');
	hold off

	% now regress out the orthogonalized phase
	c = [1 0 0 ; 0 0 0; 0 0 0];
	tmp = tmp/sum(tmp);
	X = [ tmp reg2 reg3];


	for p = (s-1)*Npix+1: s*Npix

		bhat = pinv(X)* squeeze(mag(:,p));
		mag2(:,p) = mag(:,p) - X * c *bhat;
		%plot(mag2(:,p),'r'); hold on; plot(mag(:,p)), hold off, pause
	end

end

write_nii('mag_out.nii', mag2, h,0);

h.dim(5) = 1;

vbefore = var(mag,0,1);
vafter = var(mag2,0,1);
write_nii('varBefore.nii', vbefore  , h, 0);
write_nii('varAfter.nii', vafter  , h, 0);
write_nii('rho.nii', rho*1000, h,0);

return

