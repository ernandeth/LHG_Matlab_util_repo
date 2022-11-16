function [clean junkcomps] = compcor12(dirty, hdr, Ncomps, X)
% function [clean junkcomps] = compcor12(dirty, hdr [,Ncomps] [,DesMat])
%
% written by Luis Hernandez-Garcia at UM (c) 2012
%
% the dirty input data should be a matrix : Tframes by Npixels
% The hdr is just to let us know the dimensions of the images
%
% sample usage:
%   [dirty hdr] = read_img('dirtyFiles.img');
%   [clean junkcoms] = compcor12(dirty, hdr);
%   write_img('cleanFile.img', clean, hdr);
%
% Preferred method: you can put the junkcomps into the
% design matrix of a GLM.  In that case, it may be good
% to orthogonalize the junk regressors (J) from the
% original matrix (X), like this:
%
%   for n=1:Ncomps
%       J2(:,n) = J(:,n) - X*pinv(X)*J(:,n)
%   end
%

if nargin<3
    Ncomps = 10;
end

TopPercentSigma = 30;

fprintf('\nFinding the top %d components in the noisiest pixels \n(pixels with top %d percent of the variance)\n...',Ncomps, TopPercentSigma)

if isfield(hdr,'originator')
    hdr=nii2avw_hdr(hdr);
end

sigma = std(dirty,1);
msk = ccvarmask(sigma, TopPercentSigma);
%msk = ones(size(sigma));

figure;
subplot(211)
lightbox(reshape(sigma,hdr.xdim, hdr.ydim, hdr.zdim));
title('std. dev. map')
subplot(212)
lightbox(reshape(msk,hdr.xdim, hdr.ydim, hdr.zdim));
title('pixels for compcor');
set(gcf,'Name', 'Compcor Uses the Noisiest voxels')

mdirty = dirty .* repmat(msk, hdr.tdim,1);
mdirty(isinf(mdirty))=0;
mdirty(isnan(mdirty))=0;

[u, s,v]=svd(mdirty',0);

figure
subplot(221), plot(diag(s)); title('Eigenvalues')
junkcomps = v(:,1:Ncomps);
subplot(222)
plot(junkcomps), title (sprintf('First %d components',Ncomps));
set(gcf,'Name', 'SVD identified Noise components')

% mean center the components:
junkcomps = junkcomps - repmat(mean(junkcomps,1), hdr.tdim,1);

if nargin==4
    if ~isempty(X)
        fprintf('\nExamine correlation of components with design matrix\n...')
        badinds=[];
        for n=1:size(junkcomps,2)
            % decorrelate the junk components from the design matrix (desired effects)
            % junkcomps(:,n)= junkcomps(:,n) - X*pinv(X)*junkcomps(:,n);
            
            % test the correlation between junkcomps and the design matrix.
            % if they are correlated do not use them for clean up.
            
            design =  X*pinv(X)*junkcomps(:,n);
            rho = corrcoef( junkcomps(:,n),design);
            fprintf('\nCorrelation of %d-th component with X = %f', n, rho(1,2));
            if abs(rho(1,2)) > 0.5;
                fprintf('... will not remove noise Component %d');
               badinds = [badinds ; n];
            end
            
        end
        junkcomps(:,badinds) = [];
    end
end

bhat = pinv(junkcomps)*dirty;
clean = dirty - junkcomps*bhat;
cc_sigma = std(clean,1);

% calculate the BIC for this noise model
RSS = sum(clean.^2,1);
n = size(clean,1);
k = size(junkcomps,2);
BIC = n*log(RSS/n) + k*log(n);
BIC = BIC .*msk;
BIC(BIC==0) = nan;
subplot(224)
lightbox(reshape(BIC,hdr.xdim, hdr.ydim, hdr.zdim));
title('BIC at each voxel');
subplot(223)
hist(BIC(:),50); title(sprintf('BIC histogram. Mean %f', mean(BIC(~isnan(BIC)))));

figure;
subplot(211)
lightbox(reshape(sigma,hdr.xdim, hdr.ydim, hdr.zdim));
title('std. dev. map BEFORE')
subplot(212)
lightbox(reshape(cc_sigma,hdr.xdim, hdr.ydim, hdr.zdim));
title('std. dev. map AFTER');
set(gcf,'Name', 'Compcor Reduction in Noise')


return

function msk = ccvarmask(varimage , th)
% makes a mask that preserves the data with the top TH percentage of the
% values

ordered = sort(varimage(:));
Nintgrl = cumsum(ordered)/sum(ordered(:)) * 100;
thval = find(Nintgrl>100-th);
%subplot(211), plot(ordered); title('Ordered Std. Deviations')
%subplot(212), plot(Nintgrl); title('Intrageted , Normalized Std. Deviations')

thval = ordered(thval(1))
msk = varimage;
msk(msk < thval) = 0;
msk(msk>=thval) = 1;
return
