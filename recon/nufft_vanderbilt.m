% This is code by Anuj Sharma from Vanderbilt:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate the lookup table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load spiralk;
Nk = length(k);
N = [128 128]; % image size
os = 2; % k-space oversampling factor
w = 7; b = 16.2;
for kk = 1:Nk
    mat = zeros(os*N);
    for ii = -(w-1)/2:(w-1)/2
        for jj = -(w-1)/2:(w-1)/2
            % find grid point in a w-point neighborhood
            km = round(os*k(kk,1) + ii);
            kn = round(os*k(kk,2) + jj);
            kwm = 1/w*besseli(0,b*sqrt(1-((os*k(kk,1)-km)./(w/2)).^2));
            kwn = 1/w*besseli(0,b*sqrt(1-((os*k(kk,2)-kn)./(w/2)).^2));
            kw = kwm.*kwn;
            % convert to matrix indices
            indsm = km + os*N(1)/2 + 1;
            indsn = kn + os*N(2)/2 + 1;
            if indsm > 0  && indsm <= os*N(1) && indsn > 0  && indsn <= os*N(2)
                mat(indsn,indsm) = kw;
            end
        end
    end
    table{kk} = sparse(mat);
end
save nuffttable table;

%%%%%%%%%%%%
Recon. code
%%%%%%%%%%%%

load rt_spiral_03_norm % 6-shot spiral data
% Compute KB kernel
w = 7; b = 16.2; os = 2; N = [128 128];
[wm wn] = meshgrid(-(w-1)/2:(w-1)/2);
kwm = 1/w*besseli(0,b*sqrt(1-(wm./(w/2)).^2));
kwn = 1/w*besseli(0,b*sqrt(1-(wn./(w/2)).^2));
kw = zeros(os*N);
kw(os*N(1)/2-floor(w/2)+1:os*N(1)/2+ceil(w/2),...
    os*N(2)/2-floor(w/2)+1:os*N(2)/2+ceil(w/2)) = kwm.*kwn;

p.kw = kw;
p.N = N;
p.os = os;
p.k = k;
niter = 1:5;
for ii = 1:length(niter)
    x = lsqr(@(x,transp)afun(x,transp,p),d,[],niter(ii));
    img(:,:,ii) = reshape(x,N);
end


function y = afun(x,transp,p)
if strcmp(transp,'transp')
  % call type I nufft (non-Cartesian -> Cartesian)
  y = kbnufft2d_typeI(x,p);
else
  % call type II nufft (Cartesian -> non-Cartesian)
  y = kbnufft2d_typeII(x,p);
end


function y = kbnufft2d_typeI(x,p)
% x = non-uniform data
% p = parameters structure
% OUT: y = NUFFT'd signal in uniform domain

% Load lookup-table for the 6-shot spiral
load nuffttable;

os = p.os;
N = p.N;
kw = p.kw;
Nk = length(x);
kwf = ifftshift(ifft2(ifftshift(kw)));

dC = zeros(os*N);
for ii = 1:Nk
   dC = dC + x(ii)*table{ii};
end
% Take IFFT to go to image space
y = ifftshift(ifft2(ifftshift(dC)));
% Deapodize
kwf = ifftshift(ifft2(ifftshift(kw)));
y = y./kwf;
% Truncate to original size
y = y(os*N(1)/2-N(1)/2+1:os*N(1)/2+N(1)/2,...
    os*N(2)/2-N(2)/2+1:os*N(2)/2+N(2)/2);
y = y(:); % vectorize



function x = kbnufft2d_typeII(y,p)
% y: image vector
% p: paramters struct
% x: non-uniform k-space samples

% Load lookup-table for the 6-shot spiral
load nuffttable;

os = p.os;
N = size(y);
N = ([sqrt(N(1)) sqrt(N(1))]);
Nk = length(table);
kw = p.kw;
k = p.k;
yy = reshape(y,N);

% Zero-pad
mask = zeros(os*N(1),os*N(2));
mask(os*N(1)/2-N(1)/2+1:os*N(1)/2+N(1)/2,...
    os*N(1)/2-N(1)/2+1:os*N(1)/2+N(1)/2) = yy;
yy = mask;

% Pre-emphasize
kwf = ifftshift(ifft2(ifftshift(kw)));
yy = yy./kwf;

% Take FFT to go to uniformly spaced k-space
dC = fftshift(fft2(fftshift(yy)))/(os*N(1)*os*N(2));
dCcon = conv2(dC,kw,'same');

% Resample onto non-uniform trajectory
for ii = 1:Nk    
    % Or pick the nearest Cartesian value
    km = round(os*k(ii,1));
    kn = round(os*k(ii,2));
    % convert to matrix indices
    indsm = km + os*N(1)/2 + 1;
    indsn = kn + os*N(2)/2 + 1;
    x(ii) = dCcon(indsn, indsm);
end

x = x(:);

