function myCompCor_mod(rootname, DesMatName, varargin)

% function myCompCor_mod(rootname, DesMatName, A (optional), corr_thresh_Z (optional), baseline (optional))
%
% (C) 2007-2008 Brandon W. Sur & Paul J. Chowdhry
% University of Michigan
% Last edit: June Aug 29, 2008
% report bugs to: pchowdhr@umich.edu
%
% INPUTS:
%   rootname1: name of uncorrected ANALYZE or NIFTI file
%   DesMatName: name of Design Matrix file, in ASCII (.txt) or binary format (.mat)
%           if in binary format, design matrix must be saved as variable
%           DesMat in the workspace (.mat) file
%   A: fraction of voxels to be used for tnoiseROI (in percentage)
%       (default, A = 2)
%   corr_thresh_Z:  Threshold for correlation between noise components ID'd
%       by PCA and the design matrix, expressed as a Z-score.  Noise
%       components below this threshold will be removed
%       (default, corr_thresh_Z = 13 -- !! effectively bypassing component
%       thresholding !!)
%   baseline:  1 if DesMat contains an explicit baseline as the first
%       column; 0 if there is no explicit baseline
%
% OUTPUT:
%   This program performs CompCor(Component Based Noise Correction) and
%   generates a file called 'CompCor_residuals.img'

warning off

% Optional arguments / default values
switch(size(varargin,2))
    case 0
        baseline = 0;
        A = 2;
        corr_thresh_Z = 12;
    case 3
        if ~isempty(varargin{1})
            A = varargin{1}
        else
            A = 2
        end
        if ~isempty(varargin{2})
            corr_thresh_Z = varargin{2}
        else
            corr_thresh_Z = 13
        end
        if ~isempty(varargin{3})
            baseline = varargin{3}
        else
            baseline = 0
        end
    otherwise
        error('You must enter all or none of the optional arguments.  If you wish to enter selected optional arguments, please enter an empty matrix [] for the other arguments so that they will go to their default values.');
end
     

% Import design matrix; differs depending if DesMatName is ASCII or binary
if isstr(DesMatName)
    x = DesMatName;
    if length(x) > 3
        suffix = x(end-3:end);
        switch(suffix)
            case '.txt'
                DesMat = load(DesMatName);
            case '.dat'
                DesMat = load(DesMatName);
            case '.mat'
                load DesMatName DesMat
            otherwise
                load DesMatName DesMat
        end
    end
else
    if ~isstr(DesMatName) & isa(DesMatName,'numeric')
        DesMat = DesMatName;
    end
end

if exist('DesMat','var') & isa(DesMat,'numeric')
    fprintf('Design Matrix successfully imported\n');
else
                A = 2;
    error('Design Matrix not successfully imported\n');
end

if length(rootname) > 3
        suffix = rootname(end-3:end);
        switch(suffix)
            case '.img'
                A = 2; rootname = rootname(1:end-4);
            case '.nii'
                rootname = rootname(1:end-4);
        end
end

[Y, h] = read_img(rootname);    
if isfield(h,'xdim')
    isAVW = 1;
    fprintf('Your data are in AVW format\n');
    h = avw2nii_hdr(h);
end

if length(rootname)> 3
    suffix = rootname(end-3:end);
    switch(suffix)
        case '.img'
            h = avw2nii_hdr(h);
    end
end


% 1. Detrend data
 fprintf('\nQuick detrend ...');
 Y1 = mydetrend(Y);

% 2. Temporal Noise ROI : identify the highest variance voxels
% don't consider the others for PCA

fprintf('\nIdentifying high variance pixels ');
NewData = tnoiseROI(Y1, h, DesMat, A);

% 3. Principal Component Analysis
%   a. On Anatomical Noise ROINcomponents (DEFUNCT)
%   b. On Temporal Noise ROI

fprintf('Doing the PCA ....');
Yclean = CompCor_PCA_wash_corrtest(Y, NewData, DesMat,[],corr_thresh_Z,baseline);
fprintf('\nWriting out clean data (CompCor_residuals.nii)');
% write_nii('CompCor_residuals.nii', Yclean, h, 0);

h_avw = nii2avw_hdr(h);
write_img('CompCor_residuals.img', Yclean, h_avw);

% 4. diagnostics : how much is the variance reduced?
fprintf('\ncalculating variance before and after');
hv = h;
hv.dim(5)=1;
var1 = (std(Y,0,1)).^2;
% write_nii('CompCor_varBefore.nii',var1, hv,0);

hv_avw = nii2avw_hdr(hv);
write_img('CompCor_varBefore.img', var1, hv_avw);

var2 = (std(Yclean,0,1)).^2;

% write_nii('CompCor_varAfter.nii',var2, hv,0);

write_img('CompCor_varAfter.img', var2, hv_avw);

fprintf('\nMoving output data into the tCompCor directory\n');
! mkdir tCompCor
! mv CompCor_* tCompCor

%%
function NewData = tnoiseROI(data, header, DesMat, upper_frac)

% function NewData = tnoiseROI(data, header, DesMat)
%
% (C) 2007 Brandon W. Sur
% University of Michigan
% Last edit: June 26, 2008 by Paul J. Chowdhry
% report bugs to: pchowdhr@umich.edu
%
% INPUTS:
%   data: uncorrected image data
%   header: header from read_img
%   DesMat: Design Matrix
%   Nframes: Number of time points in the data
%
% OUTPUT:
%   This program identifies voxels with high temporal standard deviation
%   and generates a file called 'tnoiseROI.nii', which contains only those voxels.
%   

[m n] = size(data);

NewData = data; % NewData is used to store noise voxels

% 1.Calculate the variance of the original time series.

tSTD = zeros(n,1);

for k = 1:n
    tSTD(k,1) = std(data(1:m,k));
end

% 2.Based on the threshold chosen by user, construct the tSTD noise ROI.

[sort_tSTD, I2] = sort(tSTD,'descend');

per = floor(upper_frac*n/100);
ind_noise_vox = I2(1:per,1);

[z3 z4] = size(ind_noise_vox);

for k = 1:n
    if k ~= ind_noise_vox(1:z3,1)
        NewData(:, k) = 0;
    end
end

mask1 = NewData(1,:);
mask1(find(mask1 ~=0)) = 1;


% 3.Perform correlation calculation to remove any tnoise voxels with high correlation
% with design matrix.

[RHO PVAL] = corr(NewData, DesMat);

low_Pvals_ind = find(abs(PVAL) < 0.2);
[p1 p2] = size(PVAL);
[r c] = ind2sub([p1 p2], low_Pvals_ind);
NewData(:,r) = 0;

mask2 = NewData(1,:);
mask2(find(mask2 ~=0)) = 1;

mn = var(NewData,0,1);
lightbox(reshape(mn,header.dim(2), header.dim(3), header.dim(4))) ; 
title('top variance map'); % pause(2)

lightbox(reshape(mask1 + mask2,header.dim(2), header.dim(3), header.dim(4)));
title('top variance mask')

% 4.Generate a NIFTI files .
write_nii('tnoiseROI.nii', NewData, header, 0);
header.dim(5)=1;
write_nii('varTnoiseROI.nii', mn, header, 0);

save tnoiseROI.mat

return

%%
function Yclean = CompCor_PCA_wash_corrtest(Y0, Yred, DesMat, Ncomponents, corr_thresh_Z, baseline)

% function CompCor_PCA_wash_corrtest(raw_data, reduced_data, DesMat, Ncomponents)
%
% (C) Paul Chowdhry 
% Last edits: 6/27/2008
% Based on CompCor_PCA_wash by Luis Hernandez-Garcia and Magnus Ulfarsson
%
% This program carries out a Principal Component Analyis on FMRI temporal noise time series data
% and removes temporal components from the raw data after determining the
% correlation of each component with the design matrix
%
% INPUTS:
%   Y0:  raw time series data
%
%   Yred:  time series of data with highest variance
%
%   DesMat: your design matrix (time x regs).  The program computes the correlation between
%   the component under scrutiny and your regressors (so you can determine
%   whether you are looking at noise or signal)
%
%   Ncomponents: [optional]  This is the maximum number of components that
%   the program will show you.  It defaults to half ot the dimensionality
%
% OUTPUT:
%   The End result is written to a file called 'residuals.nii' which contains
%   a cleaner version of your data
%
% It will display them spatial distribution (Z maps) of each temporal
% component in your data.
%
% Warning:  It [a, b] = size(DesMat);will only show you the top half of the eigenvectors.CompCor_PCA_wash_corrtest.m
% We figured that anythwarning off

%%Magnus' Code to compute PCA

[T,M] = size(Yred);
% remove mean from the data
mu = mean(Yred)';
Yred = Yred - ones(T,1)*mu';

df = T-1;

% compute the dimensionality of the matrix
[r] = laplace_pca(Yred);
fprintf('\nData Dimensionality: %d . Begin PCA decomposition', r);

% do the decomposition
% G = temporal components
% s2 = noise variance (residual)
% l = variance of principal components
% lambda = eigenvalues
% trSy = total variance in data
[G,s2,l, lambda,trSy] = nPCA(Yred,0,r);

if(r>64)
    fprintf('\n Warning.  More than 64 dimensions.  Reducing to 64');
    r=64;
end

P=Yred*G;

%load('PCA.mat')
Tcomponents = P;
Scomponents = G;

figure
imagesc(Tcomponents); colormap('gray'); title('time components to test'); drawnow
fprintf('\nBegin GLM to identify location of top half of the components ...');
contrasts = eye(size(Tcomponents,2));
X = DesMat;

BadRegs = [1:10];
Xnoise = Tcomponents(:,BadRegs);
Xnoise = [ones(length(Xnoise), 1) Xnoise];

if baseline
	DesMat_nb = DesMat(:,2:end);
	[a, b] = size(DesMat);
	[Xnoise_ind Z] = CorrMat(DesMat_nb,Xnoise,corr_thresh_Z);

else
	[a, b] = size(DesMat);
	[Xnoise_ind Z] = CorrMat(DesMat,Xnoise,corr_thresh_Z);

end

% estimate paramaters for components
betahat = pinv(Xnoise_ind)*Y0;
% remove components -
Yclean = Y0 - Xnoise_ind * betahat;

figure
imagesc(Xnoise_ind); colormap('gray'); title('Removing these components'); drawnow

save noise_regs.mat Xnoise_ind
save noise_regs.txt Xnoise_ind -ascii

save compcor_PCA.mat

return



%%
function [k,p] = laplace_pca(data, e, d, n)
% LAPLACE_PCA   Estimate latent dimensionality by Laplace approximation.
%
% k = LAPLACE_PCA([],e,d,n) returns an estimate of the latent dimensionality
% of a dataset with eigenvalues e, original dimensionality d, and size n.
% LAPLACE_PCA(data) computes (e,d,n) from the matrix data
% (data points are rows)
% [k,p] = LAPLACE_PCA(...) also returns the log-probability of each
% dimensionality, starting at 1.  k is the argmax of p.

if ~isempty(data)
    [n,d] = size(data);
    m = mean(data);
    data0 = data - repmat(m, n, 1);
    e = svd(data0,0).^2;
end
e = e(:);
% break off the eigenvalues which are write_nii( 'CompCor_residuals.nii'identically zero
% i = find(e < eps);
% e(i) = [];

kmax = min([d-1 n-2]);
%kmax = min([kmax 15]);
ks = 1:kmax;

% normalizing constant for the prior (from James)
% the factor of 2 is cancelled when we integrate over the 2^k possible modes
z = log(2) + (d-ks+1)/2*log(pi) - gammaln((d-ks+1)/2);
for i = 1:length(ks)
    k = ks(i);
    e1 = e(1:k);
    e2 = e((k+1):length(e));
    v = sum(e2)/(d-k);
    p(i) = -sum(log(e1)) - (d-k)*log(v);
    p(i) = p(i)*n/2 - sum(z(1:k)) - k/2*log(n);
    % compute logdet(H)
    lambda_hat = 1./[e1; repmat(v, length(e2), 1)];
    h = 0;
    for j1 = 1:k
        for j2 = (j1+1):length(e)
            h = h + log(lambda_hat(j2) - lambda_hat(j1)) + log(e(j1) - e(j2));
        end
        % count the zero eigenvalues all at once
        h = h + (d-length(e))*(log(1/v - lambda_hat(j1)) + log(e(j1)));
    end
    m = d*k-k*(k+1)/2;
    h = h + m*log(n);
    p(i) = p(i) + (m+k)/2*log(2*pi) - h/2;
end
p=real(p);
[pmax,i] = max(p);
k = ks(i);
return

%%
function [G,s2,l,lambda,trSy]=nPCA(Data,covar,r)
% -----------------------------------
% Usage: [log_lik]=nPCA(Data,covar,T,M)
% -----------------------------------
% Computes the nPCA log-likelihood.
% -----------------------------------
% Input:  covar=1: Data is the covariance matrix, else TxM data matrix
%         r:       The number of PCs
% -----------------------------------
% Output:  G:       Mxr orthonormal loading matrix
%          s2:      Noise variance
%          l:       Variance of the PCs
%          lambda = eigenvalues
%          trSy = total variance in data
% ----------------------------------
% Magnus Orn Ulfarsson, 2007.
% -----------------------------------
if(covar==1)
    Sy=Data;
    trSy=trace(Sy);
    [G,Lambda]=svd(Sy);
    G=G(:,1:r);
    lambda=diag(Lambda);
    s2=mean(lambda(r+1:end));
    l=lambda(1:r)-s2;
else
    Y=Data;
    [T,M]=size(Y);
    trSy=trace(Y*Y')/T;

    Y=Y-ones(T,1)*mean(Y);
    if(T>=M)
        [G,Lambda_sqrt]=svd(Y'/sqrt(T));
        G=G(:,1:r);
    else
        [P,Lambda_sqrt]=svd(Y/sqrt(T));
        G=1/sqrt(T)*Y'*P(:,1:r)*diag(diag(Lambda_sqrt(1:r,1:r)).^(-1));
    end
    lambda=diag(Lambda_sqrt).^2;
    s2=mean(lambda(r+1:end));
    l=lambda(1:r)-s2;
end

%%
function mkMIP(img, stdT)
% function mkMIP(img, stdT)
%
% Display crude maximum intensity projections of THRESHOLDED
% 3D brain maps
%
% img:  the  data
% stdT:  How many standard deviations to use for the threshold

% discard the zeros, and take the abs()
img = abs(img);
buffer = img(:);
buffer = buffer(find(buffer));

% basic stats on the image
m = mean(buffer(:));
sd = std(buffer(:));

threshold = m + stdT *sd;

fprintf('\nMIP info: -- Mean:  %6.2f  | StdDev: %6.2f  | Threshold: %6.2f | Max:  %6.2f | Min: %6.2f' , ...
    m, sd, threshold, max(img(:)), min(img(:)) );

% we use 10% of the mean to create a background for overlaying the
% datarmReg
buffer = zeros(size(img));
buffer(find(img)) = m/10;

img(find(img < threshold)) = 0;
img = img + buffer;

% do the MIP here
f1 = squeeze(sum(img,1));
f2 = squeeze(sum(img,2));
f3 = squeeze(sum(img,3));

colormap(jet)

set(gcf, 'Position',[100 100 800 400])
subplot (231), imagesc(f1'), axis xy, axis equal, axis tight
subplot (232), imagesc(f2'), axis xy, axis equal, axis tight
subplot (233), imagesc(f3'), axis xy, axis equal, axis tight


return

%%
function CompCor_rmReg (root, allDesMat, xclude)

% rmReg (file_root, allslices_Desmats. ,[xclude])
%
% removes known trends in time series of images.  Aimed at Physio correction
%
%   (c) 2006 Luis Hernandez-Garcia
%   University of Michigan
%   report bugs to:  hernan@umich.edu
%   Modified by Brandon Sur at UM on 06/20/2007
%
% The program solves the linear model specified in DesMat for the
% parameters by ordinary least squares. It returns the residuals to a 4D image file
% called residuals.img
%
% important note.  we estimate the intercept regressor, but we do not remove it from
% the data.  The program assumes the intercept is the FIRST column in the
% input (allDesMat).
%
% xclude lets you estimate, but not remove regressors in the
% design matrix.  If ommited we default to 1 (not removing the first
% regressor, which is assumed to be the baseline)

if nargin==2
    xclude=1;
end

[all_data, h] = read_img(root);
Nframes = size(all_data,1);

% check to see if these are AVW or NIFTI files:
if length(root)> 3
    suffix = root(end-3:end);
    switch(suffix)
        case '.img'
            h = avw2nii_hdr(h);
    end
end

output = zeros(size(all_data));
fprintf('\n Done reading the data. crunching ...');
%hnames = dir(sprintf('%s*.hdr',root));
%h = read_hdr(hnames(1).name);
Spix = h.dim(2) *h.dim(3);

slice = 0;
pcount = inf;

hv = h;
hv.dim(1) = 3;
hv.dim(5) = 1;
% compute the variance after removing regressors
var1 = (std(all_data,0,1)).^2;
write_nii('CompCor_varBefore.nii',var1, hv,0);

for pix =1:size(all_data,2)

    if pcount > Spix
        pcount=1;
        slice = slice +1;
        % check to see whether we have a different design matrix for each
        % slice...
        if size(size(allDesMat),2) == 3
            DesMat = squeeze(allDesMat(slice,:,:));
        else
            DesMat = allDesMat;NewData(:,ind_noise_vox(k,1)), DesMat(:,k)
        end
        % adjust the design matrix' size by taking only the end
        % i.e. - this clips off the physio data collected during disdaqs
        DesMat = DesMat((end-Nframes+1):end,:);
        xtx_inv = pinv(DesMat);
        fprintf('\r  slice %d ...',slice);
    end

    pix_data = all_data(:,  pix );
    beta_est = xtx_inv*pix_data;

    %    if pix==5000
    %        keyboard
    %    end
write_nii
    xDM = DesMat(:,xclude+1:end);
    xBeta = beta_est(xclude+1:end,:);
    output(:, pix ) = pix_data - xDM *xBeta ;
    pcount = pcount +1;

end

% compute the variance after removing regressors
var2 = (std(output,0,1)).^2 ;
write_nii('CompCor_varAfter.nii',var2,hv,0);

fprintf('\n Writing output residuals ...');
outh = h;
outh.dim(5) = size(all_data,1);


%write_hdr( 'residuals.hdr', outh);
%write_img( 'residuals.img', output, outh);
write_nii( 'CompCor_residuals.nii', output, outh,0);

warning on

save regressors xDM DesMat allDesMat

fprintf('\n...Done');
return

%%
function result = mydetrend( data, saveCoeffs)
%function result = mydetrend( data [,saveCoeffs])
% this function fits and subtracts a third order polynomial
% from the data.
%
% (c) 2005 Luis Hernandez-Garcia
% University of Michigan
% report bugs to:  hernan@umich.edu
%
% Modified by Brandon Sur on July 5, 2007
%
%
[m n] = size(data);
t = linspace(0, m, m);
t =t';

for k = 1:n
    [coeffs, error] = polyfit(t,data(:,k),3);
    data(:,k) = data(:,k)-polyval(coeffs,t);
end

result = data;

if nargin<2
    save coeffs
end
return

%%
function spmJr (root, DesMat, contrast)

% spmJr (file_root, Desmat, contrast_matrix)
%
% The extremely quick and dirty SPM analysis
%
%   (c) 2006 Luis Hernandez-Garcia
%   University of Michigan
%   report bugs to:  hernan@umich.edu
%
% The program solves the linear model specified in DesMat for the
% parameters by ordinary least squares.  It calculates T scores and Z scores
%
% I'm borrowing (calling) spm_t2z.m  for the Z score calculation.  Thanks!
%
% No filtering, whitening, or
% anything at all is done to the data.  You're on your own for that.
%
% inputs:
%
% file_root:  this is where the data live:  either a single 4D AVW format or
%             a time series of files containing the individual frames.
%             if this is not a string, I assume you enetered a t x p matrix, where
%             columns are pixels and rows are time points
% DesMat:     The design matrix.  The dimensions must be [time x regressors]
% contrast:   a matrix of contrasts.  Each contrast is a row of the matrix.
%
%
% outputs:
%
% Tmap_????.img   : The t-test for each contrast specified in the matrix
%                  (rows) evaluated at each voxel.  It's in AVW format.
% Zmap_????.img   : the corresponding Z score associated with the t test
% spmJr.mat       : A .mat file containing the following variables used in
%                   the analysis:
%               - beta_est (the parameter estimates for each regressor)
%               - DesMat (the design matrix used in the estimation)
%               - contrast (the contrast matrix)
%               - df (the degrees of freedom in the Z score calculation)
%               - vCon (the estimated variance of the contrast of parameter estimates.
% That sounds more complicated than it is, doesn't it?  You've estimated a parameter,
% ie - the amplitude corresponding to a regressor.  This is a measure of
% the variance or uncertainty of that estimate and it's calculated from the
% variance of the data at that voxel and from the design matrix
%

names = dir(sprintf('%s*',root));
if isstr(root)
    names = dir(sprintf('%s*',root));
    suffix = names(1).name(end-3:end);
    if suffix=='.nii'
        [all_data, h] = read_img(names(1).name);
        h = nii2avw_hdr(h);
    else
        all_data = read_img(root);
        hnames = dir(sprintf('%s*.hdr',root))
        h = read_hdr(hnames(1).name);
    end
    fprintf('\n Done reading the data. crunching ...');
else
    all_data=root;
end

Ncon = size(contrast,1);

tmap = zeros( Ncon, size(all_data, 2));
vCon = zeros( Ncon, size(all_data, 2));
Tol = 10^(-16);
df = size(all_data,1) - size(DesMat,2)+1;

warning off

fprintf('\n Estimating Beta parameters and variance ...');

xtx_inv = pinv(DesMat);
beta_est = xtx_inv*all_data;
V = var(all_data);

for n=1:size(contrast,1)
    fprintf('\n Contrast n. %d ...',n);

    for pix=1:size(all_data,2)
        vCon(n,pix) = contrast(n,:) * xtx_inv * V(pix) * xtx_inv' * contrast(n,:)';

    end

    tmap(n,:) = (beta_est' * contrast(n,:)') ./ sqrt(vCon(n,:))';
end

tmap(find(isnan(tmap)))=0;
tmap(find(isinf(tmap)))=0;
zmap = tmap;
zmap = spm_t2z(tmap(:),df);
zmap = reshape(zmap,Ncon, size(all_data, 2));

fprintf('\n Writing output files (Tmap, Zmap) ...');

if isstr(root)
    for n=1:Ncon
        % make sure that we write stats maps as floats.
        outh = h;
        outh.tdim = 1;
        outh.datatype=16;
        outh.bits = 32;
        outh.glmax = max(tmap(n,:));
        outh.glmin = min(tmap(n,:));

        write_hdr( sprintf('Tmap_%04d.hdr', n), outh);
        write_img_data( sprintf('Tmap_%04d.img',n), tmap(n,:), outh);

        outh.glmax = max(vCon(n,:));
        outh.glmin = min(vCon(n,:));
        write_hdr( sprintf('ConVar_%04d.hdr', n), outh);
        write_img_data( sprintf('ConVar_%04d.img',n), vCon(n,:), outh);

        outh.glmax = max(zmap(n,:));
        outh.glmin = min(zmap(n,:));
        write_hdr( sprintf('Zmap_%04d.hdr',n), outh);
        write_img_data( sprintf('Zmap_%04d.img',n), zmap(n,:), outh);
        warning on
    end
end

fprintf('\nSaving Intermediate data to spmJr.mat file ...');
save spmJr beta_est contrast DesMat vCon df zmap tmap

fprintf('\n...Done');
return

%%
function hdr = nii2avw_hdr(niih)
%funtion avwh = nii2avw_hdr(niih)

hdr = define_avw_hdr;

hdr.sizeof_hdr = niih.sizeof_hdr;
%hdr.pad1 = niih.?
hdr.extents= niih.extents;
%hdr.pad2  = niih.?
%hdr.regular= niih.?
%hdr.pad3  = niih.?
hdr.dims = niih.dim(1);
hdr.xdim = niih.dim(2);
hdr.ydim = niih.dim(3);
hdr.zdim = niih.dim(4);
hdr.tdim = niih.dim(5);
%hdr.pad4 = niih.?
hdr.datatype = niih.datatype;
hdr.bits  = niih.bitpix;
%hdr.pad5 = niih.?
hdr.xsize = niih.pixdim(2);
hdr.ysize = niih.pixdim(3);
hdr.zsize = niih.pixdim(4);
%hdr.pad6 = niih.?
hdr.glmax = niih.cal_max;
hdr.glmin = niih.cal_min;
hdr.descrip = niih.descrip;
hdr.aux_file = niih.aux_file;
%hdr.orient = niih.?
%hdr.origin = ?
%hdr.generated = niih.?
%hdr.scannum = niih.?
%hdr.patient_id = niih.?
%hdr.exp_date = niih.?
%hdr.exp_time = niih.?
%hdr.hist_un0 = niih.?
%hdr.views = niih.?
%hdr.vols_added = niih.?
%hdr.start_field = niih.?
%hdr.field_skip = niih.?
%hdr.omax = niih.?
%hdr.omin = niih.?
%hdr.smax = niih.?
%hdr.smin = niih.?

return

%%
function hdr = define_avw_hdr()
% function define_avw_hdr()
% creates a structure with a blank analyze header'
hdr = struct(...
    'sizeof_hdr'      , 0, ...
    'pad1'            , char(zeros(28,1)), ...
    'extents'         , 0 , ...
    'pad2'            , char(zeros(2,1)), ...
    'regular'         , 'r', ...
    'pad3'            , ' ',...
    'dims'            , 0, ...
    'xdim'            , 0, ...
    'ydim'            , 0, ...
    'zdim'            , 0, ...
    'tdim'            , 0, ...
    'pad4'            , char(zeros(20,1)),...
    'datatype'        , 0, ...
    'bits'            , 0, ...
    'pad5'            , char(zeros(6,1)),...
    'xsize'           , 0, ...
    'ysize'           , 0, ...
    'zsize'           , 0, ...
    'pad6'            , char(zeros(48,1)) ,...
    'glmax'           , 0, ...
    'glmin'           , 0, ...
    'descrip'         , char(zeros(80,1)),...
    'aux_file'        , char(zeros(24,1)),...
    'orient'          , char(zeros(1,1)), ...
    'origin'          , zeros(5,1),...
    'generated'       , char(zeros(10,1)),...
    'scannum'         , char(zeros(10,1)),...
    'patient_id'      , char(zeros(10,1)),...
    'exp_date'        , char(zeros(10,1)),...
    'exp_time'        , char(zeros(10,1)),...
    'hist_un0'        , char(zeros(3,1)),...
    'views'           , 0, ...
    'vols_added'      , 0, ...
    'start_field'     , 0, ...
    'field_skip'      , 0, ...
    'omax'            , 0, ...
    'omin'            , 0, ...
    'smax'            , 0,...
    'smin'            , 0 ...
    );

return

%%

function [Xnoise_ind Z] = CorrMat(DesMat, Xnoise, Zthresh)

% function [Xnoise2 Z] = CorrMat(DesMat, Xnoise, Zthresh)

% (C) 2008 Paul J. Chowdhry
% University of Michigan
% Last edit: Aug 5, 2008
% report bugs to: pchowdhr@umich.edu
%
% INPUTS:
%   DesMat: Design Matrix
%   Xnoise: noise matrix
%   Zthresh:  Z-score threshold, below which components from Xnoise will
%       become part of Xnoise2
% OUTPUT:
%   Xnoise2: contains only columns of Xnoise whose correlation with the
%       design matrix (DesMat), measured as a Z-score, is below Zthresh
%   Z:  vector of Z-scores corresponding to the correlation of each columnx
%       of Xnoise with the design matrix (DesMat)

N = size(DesMat,1);
p = size(DesMat,2);

contrast = ones(1,size(DesMat,2));

% Code below from spmJr.m

all_data = Xnoise;

Ncon = size(contrast,1);

tmap = zeros( Ncon, size(all_data, 2));
vCon = zeros( Ncon, size(all_data, 2));
Tol = 10^(-16);
df = size(all_data,1) - size(DesMat,2)+1;

warning off

fprintf('\n Estimating Beta parameters and variance ...');

xtx_inv = pinv(DesMat);
beta_est = xtx_inv*all_data;
N = size(all_data,1);
p = size(beta_est,1);
V = var(all_data);

for n=1:size(contrast,1)
    fprintf('\n Contrast n. %d ...',n);

    for pix=1:size(all_data,2)
        %vCon(n,pix) = contrast(n,:) * xtx_inv * V(pix) * xtx_inv' * contrast(n,:)';
		var_est = (all_data(:,pix) - DesMat*beta_est(:,pix))'*...
			(all_data(:,pix) - DesMat*beta_est(:,pix)) / (N-p);
%		vCon(n, pix) = contrast(n,:) * xtx_inv * var_est * xtx_inv' * contrast(n,:)';
 		vCon(n, pix) = contrast(n,:) * inv(DesMat'*DesMat)*contrast(n,:)'*var_est;
    end
    
    tmap(n,:) = (beta_est' * contrast(n,:)') ./ sqrt(vCon(n,:))';
end

t = tmap;

df = size(Xnoise,1) - size(DesMat,2) - 1; % degrees of freedom for t2z

[Z a] = spm_t2z(t,df);

badregs = find(abs(Z)<Zthresh);

% noise matrix of components independent of the design matrix 
Xnoise_ind = Xnoise(:,badregs);

return