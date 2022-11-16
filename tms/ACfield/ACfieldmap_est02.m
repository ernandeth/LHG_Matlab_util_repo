%function ACfield_est = ACfieldmap_est(rawsignal, M0, Boff, fovxy, freq, delay, aqbw)

xdim=32;

M0 = phantom(xdim);
[xdim ydim] = size(M0);

freq = 1e3; % ACfield's frequency
aqbw = 12e3;  % acquisition bandwith in Hz
rawsignal = fft2(M0);  % observed signal
fovxy = 20; % field of view


delay = 1e-4;

[xdim ydim] = size(M0);
[kxdim kydim] = size(rawsignal);

gambar = 42.57e3;        % gamma/2pi in kHz/T
gam = gambar*2*pi;

% Generate space sampling location grid:
xmax = fovxy/2;
ymax = fovxy/2;
xstep = fovxy/size(M0,2);
ystep = fovxy/size(M0,1);
[x y] = meshgrid(-xmax:xstep: xmax-xstep  ,  -ymax:ystep: ymax-ystep );
x = x(:);
y = y(:);

% Generate a K-space sampling grid:
deltax = fovxy/xdim;
kmax  = pi/deltax;
kstep = 2*pi/fovxy;
[kx ky] = meshgrid(-kmax:kstep: kmax-kstep  ,  -kmax:kstep: kmax-kstep );
kx = kx(:)';
ky = ky(:)';

% Generate a spatial magnetic field pattern
Bcoeffs = ...
   [0 1e-6  0 0 1e-7; 
    0 0 0 0 1e-7; 
    0 0 0 0 0; 
    0 1e-6 0 0 1e-7];

ACfield = zeros(size(M0(:)));
Ncoeffs = length(Bcoeffs);
for n=1:Ncoeffs
    ACfield = ...
        ACfield + ...
        Bcoeffs(1,n)*cos((n-1)*x/pi) + ...
        Bcoeffs(2,n)*sin((n-1)*x/pi) + ...
        Bcoeffs(3,n)*cos((n-1)*y/pi) +...
        Bcoeffs(4,n)*sin((n-1)*y/pi);
end

BoffAC = ACfield;

% A field that is NOT defined by a basis set:
BoffAC = 1e-5*exp(-((x-2).^2 + (y-2).^2)/5);

imagesc(reshape(BoffAC,xdim,ydim)); 


M0 = M0(:);
threshold = 0.02* mean(M0);
inds = find(abs(M0)>=threshold);

rawsignal = rawsignal(:).';

t = linspace(0, kxdim/aqbw, kxdim);
t = repmat(t, 1, kydim);

M0phase = zeros(length(M0), length(rawsignal), 'single');
% start building  the parts of the signal equation integral:
% this is part with spin density times the the desired k-space trajectories:
% AND any DC offset in the field

F = zeros(length(M0(:)), length(rawsignal(:)), 'single');
Fi = zeros(length(M0(:)), length(rawsignal(:)), 'single');
for r=1:length(x)
    F(r,:) = exp(-i * (kx * x(r) + ky * y(r)) ) ;% .* exp (i * gambar * BoffDC(r) * t);
    Fi(r,:) = exp(i * (kx * x(r) + ky * y(r)) ) ;
end



% the next part of the equation is the phase gained because of the AC field
% this phase will be weighted by the field's amplitude at each pixel
% when we update the AC field map
ACfun = sin(linspace(-2*pi,2*pi,32));

ACfun2 = repmat(ACfun, kydim,1);
ACphase = gambar * cumsum( ACfun2,2 );
ACphase = ACphase(:).';

F_ac = exp(i * kron(BoffAC(:), ACphase));

% Next, the dot-multiplication makes sure that each pixel's k-space trajectory includes the
% external field distorsion
F2 = (F .* F_ac).';
rawsignal = F2 * M0(:);

im_distorted = reshape(Fi*rawsignal, xdim,ydim) ;
im_distorted = fftshift(ifft2(fftshift(reshape(rawsignal,xdim,ydim))));


%%
F2_est = pinv(M0(:)')*rawsignal';
F2_est = lsqr(M0, rawsignal);
%%

figure(1), subplot(211)
plot(ACfun)
title('time function')
subplot(212)
imagesc(reshape(BoffAC, xdim,ydim))
title('spatial pattern')

figure(2)
imagesc(angle(F));
title('FT matrix')


figure(3)
imagesc(angle(F_ac));
title('distorted FT matrix')

figure(4), subplot(211)
imagesc(abs(reshape(rawsignal,xdim,xdim)));
title('raw signal')

figure(4), subplot(212)
imagesc(abs(im_distorted))
title('distorted, reconned image')

 
close all;


parms.M0 = reshape(M0,xdim,ydim);
parms.ACfun = ACfun;
parms.F = F;
parms.fovxy = fovxy;

dummy = [];

tstart = tic;

Bcoeffs0 = zeros(4,10);
Bcoeffs0 = zeros(2,10);


optvar=optimset('lsqnonlin');
optvar.Display = 'iter';
% optvar.MaxIter = 50;
% optvar.MaxFunEvals = (length(rawsignal(:))+1)*400;
% optvar.TolFun = 1e-10;


%{
tstart = tic;
% call the fitting routine:
[Bcoeffs_est res] = lsqnonlin(@ACfield_signal_lsq02, Bcoeffs0,...
    [],[],...
    optvar,...
    parms,...
    rawsignal(:));

fprintf('\nTime Elapsed = %f',toc(tstart));

ACfield = zeros(size(M0(:)));
% Ncoeffs = length(Bcoeffs_est);
% for n=1:Ncoeffs
%     ACfield = ...
%         ACfield + ...
%         Bcoeffs_est(1,n)*cos((n-1)*x/pi) + ...
%         Bcoeffs_est(2,n)*sin((n-1)*x/pi) + ...
%         Bcoeffs_est(3,n)*cos((n-1)*y/pi) +...
%         Bcoeffs_est(4,n)*sin((n-1)*y/pi);
% end

% alternative code using polynomial bases
for n=1:Ncoeffs
    ACfield = ...
        ACfield + ...
        Bcoeffs_est(1,n)* (x/xmax).^(n-1) + ...
        Bcoeffs_est(2,n)* (y/ymax).^(n-1) ; 
end
subplot(211)
imagesc(reshape(ACfield,xdim,ydim)); 
title('Estimate')
subplot(212)
imagesc(reshape(BoffAC,xdim,ydim)); 
title('original')
%}

