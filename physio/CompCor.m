function CompCor(rootname, DesMat, TR, A)

% function CompCor(rootname1, DesMat, TR, A)
%
% (C) 2007 Brandon W. Sur
% University of Michigan
% Last edit: July, 12 2007
% report bugs to: wsur@umich.edu
%
% INPUTS:
%   rootname1: uncorrected ANALYZE or NIFTI file 
%   DesMat: Design Matrix
%   Nframes: Number of time points in the data
%   Nslices: Number of slices obtained at one time point
%   TR: TR is TR
%   A: fraction of voxels to be used for tnoiseROI(in percentage)
%
% OUTPUT:
%   This program performs CompCor(Component Based Noise Correction) and
%   generates a file called either 'anoiseROI_residuals.nii' or
%   'tnoiseROI_residuals.nii'.


[Y, h] = read_img(rootname);
Y1 = mydetrend(Y);

% 1. Check if the root file is an ANALYZE or NIFTI file

if length(rootname)> 3
   suffix = rootname(end-3:end);
   switch(suffix)
   case '.img'
       h = avw2nii_hdr(h);
   end
end

% 2. Anatomical Noise ROI

% anoiseROI(data, header, rawDM);

% 3. Temporal Noise ROI

tnoiseROI(Y1, h, DesMat, A);

% 4. Principal Component Analysis
%   a. On Anatomical Noise ROI
%   b. On Temporal Noise ROI

% CompCor_PCA_wash(rootname1,'anoiseROI.nii', DesMat, TR, []);
% mkdir aCompCor
% !mv CompCor_residuals.nii anoiseROI_residuals.nii
% !mv anoiseROI_residuals.nii aCompCor
% !mv CompCor_* aCompCor

CompCor_PCA_wash(rootname,'tnoiseROI.nii', DesMat, TR,[]);
! mkdir tCompCor
! mv CompCor_* tCompCor

%%
function tnoiseROI(data, header, DesMat, upper_frac)

% function tnoiseROI(data, header, DesMat)
%
% (C) 2007 Brandon W. Sur
% University of Michigan
% Last edit: July 12, 2007
% report bugs to: wsur@umich.edu
%
% INPUTS:
%   data: uncorrected image data  
%   header: header from read_img
%   DesMat: Design Matrix
%   TR: TR is TR
%   Nframes: Number of time points in the data
% 
% OUTPUT:
%   This program identifies voxels with high temporal standard deviation 
%   and generates a file called 'tnoiseROI.nii', which contains only those voxels.


[m n] = size(data);

NewData = data; % NewData is used to store noise voxels

% 1.Calculate the variance of the original time series.

tSTD = zeros(n,1);

for k = 1:n
    tSTD(k,1) = std(data(1:m,k));
end

% 2.Based on the threshold chosen by user, construct the tSTD noise ROI. 

% upper_frac = input('Enter the number for the fraction of voxels to be used for tnoiseROI(in percentage): ');

[sort_tSTD, I2] = sort(tSTD,'descend');

per = floor(upper_frac*n/100);
ind_noise_vox = I2(1:per,1);

[z3 z4] = size(ind_noise_vox);

for k = 1:n
    if k ~= ind_noise_vox(1:z3,1)
        NewData(:, k) = 0;
    end
end

% 3.Perform correlation calculation to remove any tnoise voxels with high correlation
% with design matrix.

[a, b] = size(DesMat); 

for k = 1:b
        [RHO PVAL] = corr(NewData(:,ind_noise_vox(k,1)), DesMat(:,k));
        if PVAL < 0.2 % exlcude voxels with a p-value less than 0.2
            NewData(:,ind_noise_vox(k,1)) = 0;
        end 
end

% 4.Generate a NIFTI file.

write_nii('tnoiseROI.nii', NewData, header, 0);


%%
function CompCor_PCA_wash(raw_root, rootname, DesMat, TR , Ncomponents)

% function CompCor_PCA_wash(raw_root, rootname, DesMat, TR [, Ncomponents])
% 
% (C) Luis Hernandez-Garcia and Magnus Ulfarsson at UM
% Last edits: 2/28/2007
% Modified by: Brandon Sur at UM 06/20/2007
%
% This program carries out a Principal Component Analyis on FMRI anatomical or temporal noise time series data
% and removes user-specified temporal components from the raw data
%
% INPUTS:
%   rootname:  root name of time time series data:  e.g. 'ravol'
%
%   DesMat: your design matrix (time x regs).  The program computes the correlation between
%   the component under scrutiny and your regressors (so you can determine
%   whether you are looking at noise or signal)
%
%   TR:  the sampling period of the experiment (TR)
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
% Warning:  It will only show you the top half of the eigenvectors.  
% We figured that anything you want to remove will be in there.  
% It's easy to edit that code, though.
%
% The maps are displayed as maximum intensity
% projections (think Glass Brains in SPM)
% 
% PCA_wash prints out some important information about each component:
% 1. some info about the spatial maps and how the MIPs are made
% 2. how correlated each regressor is with your design matrix
% 3. How much variance each component accounts
%
% The program will then ask you which components you want to remove.
% It will then regress them out of the data
%
% You don't want to remove stuff that is correlated with the design matrix
% 

warning off

[Y, h] = read_img(rootname);

%Magnus' Code to compute PCA
[T,M]=size(Y);

% remove mean from the data
mu=mean(Y)';
Y=Y-ones(T,1)*mu';

df = T-1;

% compute the dimensionality of the matrix
[r] = laplace_pca(Y);
fprintf('\nData Dimensionality: %d . Begin PCA decomposition', r);

% do the decomposition
% G = temporal components
% s2 = noise variance (residual)
% l = variance of principal components
% lambda = eigenvalues
% trSy = total variance in data
[G,s2,l, lambda,trSy] = nPCA(Y,0,r);

if(r>64)
    fprintf('\n Warning.  More than 64 dimensions.  Reducing to 64');
    r=64; 
end

P=Y*G;
save PCA.mat

%load('PCA.mat')
Tcomponents = P;
Scomponents = G;

figure
imagesc(Tcomponents); colormap('gray'); title('time components to test'); drawnow

fprintf('\nBegin GLM to identify location of top half of the components ...');
contrasts = eye(size(Tcomponents,2));

% Decide on the number of components to show the user
% The default is half of the dimensionality
Ncontrasts = round(size(contrasts,1)/2);
if nargin==4
    if Ncontrasts > Ncomponents
        Ncontrasts = Ncomponents;
    end
end

% Test the General Linear Model
spmJr(rootname, Tcomponents, contrasts(1:Ncontrasts,:));

maps = dir('Zmap*.img');
if ~isempty(DesMat)
    %DesMat=DesMat(:,2);
    DesMat = DesMat - ones(size(DesMat,1),1)*mean(DesMat);
end
t = 0:(size(Tcomponents,1))-1;
t = t*TR;
omega = linspace(-1/(2*TR), 1/(2*TR), size(Tcomponents,1))
for f=1:length(maps)
    %map = squeeze(G(:,f));
    %m = reshape(map, h.xdim, h.ydim, h.zdim);

    m = read_img2(maps(f).name);
    
    comp = Tcomponents(:,f);
    % compute how correlated these "Bad regressors" are with the design
    % matrix - making sure we don't throw out the baby with the bath water.
    if ~isempty(DesMat)
        correlations = abs((comp'* DesMat) ./ sqrt( trace(DesMat'*DesMat) * (comp'*comp)) );

        fprintf('\nCorrelation to the Design Matrix: %s', num2str(correlations));
        fprintf('\nCorresponding Eigenvalue (percent of the variace): %7.2f ', 100*lambda(f)/trSy);
    end

    figure
    fprintf('\n--\ncomponent N. %d', f);
    mkMIP(m, 3.5);
    subplot(223), plot(t, Tcomponents(:,f)); axis tight
    
    title(['Time Component N.  ' num2str(f)])
    subplot(224), plot(omega,abs(fftshift(fft(Tcomponents(:,f)))));grid;axis tight
    title(['FFT of Component N.  ' num2str(f)])
    drawnow
    pause
end

BadRegs = input('Enter the component numbers you wish to remove (as a vector[]):  ');

removeThis = Tcomponents(:, BadRegs);
removeThis = [ones(length(removeThis), 1) removeThis];
close all
CompCor_rmReg(raw_root, removeThis);

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
% break off the eigenvalues which are identically zero
i = find(e < eps);
e(i) = [];

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

% we use 10% of the mean to create a background for overlaying the data
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
            DesMat = allDesMat;
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








    
