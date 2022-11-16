function PCA_wash(rootname, DesMat, TR , Ncomponents)
% function PCA_wash(rootname, DesMat, TR [, Ncomponents])
% 
% (C) Luis Hernandez-Garcia and Magnus Ulfarsson at UM
% Last edits: 2/28/2007
%
% This program carries out a Principal Component Analyis on FMRI time series data
% and removes user-specified temporal components from the data
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

[Y, h] = read_img_series(rootname);

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
djsspmJr(rootname, Tcomponents, contrasts(1:Ncontrasts,:));

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
    subplot(224), plot(omega,abs(fftshift(fft(Tcomponents(:,f)))));axis tight
    title(['FFT of Component N.  ' num2str(f)])
    drawnow
    pause
end

BadRegs = input('Enter the component numbers you wish to remove (as a vector[]):  ');

removeThis = Tcomponents(:, BadRegs);
removeThis = [ones(length(removeThis), 1) removeThis];
close all
rmReg(rootname, removeThis);

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
    