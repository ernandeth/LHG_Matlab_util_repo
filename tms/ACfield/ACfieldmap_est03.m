%function ACfield_est = ACfieldmap_est(rawsignal, M0, Boff, fovxy, freq, delay, aqbw)

xdim=32;

M0 = phantom(xdim);
[xdim ydim] = size(M0);

freq = 1e3; % ACfield's frequency in Hz
aqbw = 125e3;  % acquisition bandwith in Hz
rawsignal = fft2(M0);  % observed signal
fovxy = 20; % field of view


delay = 1e-4;

[xdim ydim] = size(M0);
[kxdim kydim] = size(rawsignal);

gambar = 42.57e6;        % gamma/2pi in Hz/T
gam = gambar*2*pi;       % gamma in (rad/s) / T

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
% Ncoeffs = length(Bcoeffs);
% for n=1:Ncoeffs
%     ACfield = ...
%         ACfield + ...
%         Bcoeffs(1,n)*cos((n-1)*x/pi) + ...
%         Bcoeffs(2,n)*sin((n-1)*x/pi) + ...
%         Bcoeffs(3,n)*cos((n-1)*y/pi) +...
%         Bcoeffs(4,n)*sin((n-1)*y/pi);
% end

BoffAC = ACfield;

% A field that is NOT defined by a basis set:
BoffAC = 1e-4*exp(-((x).^2 + (y).^2)/100);

% test with uniform field:
% BoffAC(:)=1e-8;

imagesc(reshape(BoffAC,xdim,ydim)); 


M0 = M0(:);
threshold = 0.02* mean(M0);
inds = find(abs(M0)>=threshold);

dt = 1/aqbw;


M0phase = zeros(length(M0), length(rawsignal));
% start building  the parts of the signal equation integral:
% this is part with spin density times the the desired k-space trajectories:
% AND any DC offset in the field

F = zeros(length(M0(:)), length(rawsignal(:)));
Fi = zeros(length(M0(:)), length(rawsignal(:)));
for r=1:length(x)
    F(r,:) = exp(-j * (kx * x(r) + ky * y(r)) ) ;% .* exp (i * gambar * BoffDC(r) * t);
    Fi(r,:) = exp(j * (kx * x(r) + ky * y(r)) ) ;
end



% the next part of the equation is the phase gained because of the AC field
% this phase will be weighted by the field's amplitude at each pixel
% when we update the AC field map
ACfun = sin(linspace(0,0.1*pi,32));

% Test with a DC field :
%  ACfun(:) = 1e-3;

ACfun2 = repmat(ACfun, kydim,1);
ACphase = gam * cumsum( ACfun2,2 )* dt ;

ACphase = ACphase';
ACphase = ACphase(:)';

% F_ac = exp(i * kron(BoffAC(:), ACphase));
% This is equivalent to:
F_ac = exp(-j * (BoffAC(:) * ACphase));

% Next, the dot-multiplication makes sure that each pixel's k-space trajectory includes the
% external field distorsion
F2 = (F .* F_ac).';
rawsignal = F2 * M0(:);
rawsignal_clean =  F*M0(:);


im_distorted = fftshift(ifft2(fftshift(reshape(rawsignal,xdim,ydim))));
im_clean = fftshift(ifft2(fftshift(reshape(rawsignal_clean,xdim,ydim))));

figure(1), subplot(311)
plot(ACfun)
title('time function')


subplot(312)
imagesc(reshape(BoffAC, xdim,ydim))
title('spatial pattern')


figure(2)
imagesc(angle(F));
title('FT matrix')


figure(3)
imagesc(angle(F_ac));
title('distorsion added to FT matrix')

figure(4), subplot(211)
imagesc(abs(reshape(rawsignal,xdim,xdim)));
title('raw signal')

figure(4), subplot(212)
imagesc(abs(im_distorted ))
title('distorted, reconned image')

drawnow
%{
F2_est = pinv(M0(:)')*rawsignal';
F2_est = lsqr(M0, rawsignal);
%}

% This is for the DC case, where the phase difference gives you the map.
Bdc = -angle( im_distorted ./im_clean);
figure(5)
mask = phantom(xdim); mask(mask>0)=1;
imagesc(Bdc .* mask)
title('phase difference')
drawnow

% The estimation is done here:
parms.Fn = F(1,:);
parms.ACfun = ACphase;
parms.M0 = M0;
   
Br_est = zeros(size(BoffAC));
for n=1:length(BoffAC(:))   
   fprintf('\nEvaluating ... %d  BoffAC(n) = %f  ', n, 1e4*BoffAC(n));
   parms.Fn = F(:,n);
   
%   Br_est(n) =  fzero(@(B) ACfield_signal_lsq03(B, parms, rawsignal ),0 );
   
    % Not trying to estimate - check if the true value satisfies the
    % equation

    % for testing only!!
    
    % yn=sum( F(:,n) .* exp(-j .* BoffAC(:) .* ACphase(n)) .* M0(:));
    %
    % Br_est(n) =  ACfield_signal_lsq03(BoffAC(n), parms, yn );
    % fprintf('LHS-RHS = %f  ', 1e4*Br_est(n));
      
    
    % Br_est(n) =  ACfield_signal_lsq03(BoffAC(n), parms, rawsignal(n) );  
    % fprintf(' LHS - RHS = %f  ', Br_est(n));
    

    num=0;
    den=0;
    for k=1:length(M0)
        tmp1 = 0;
        tmp2 = 0;
        for kk=1: length(M0)
            tmp1 = tmp1 + F(n,k)*F(n,kk)*ACphase(k)*M0(k)*M0(kk);
            tmp2 = tmp2 + F(n,k)*F(n,kk)*(ACphase(k)^3)*M0(k)*M0(kk);
        end
        num = num + rawsignal(n)*F(n,k).*ACphase(k) - (1/2)*tmp1;
        den = den + rawsignal(n)*F(n,k).*(ACphase(k))^3 - (1/12)*tmp2;
    end
    Br_est(n) = sqrt( num / den);
    fprintf(' Boff_AC_estimate(n) : %f ', Br_est(n));
end



figure(1) ; subplot(313)
imagesc(reshape(abs(Br_est), xdim,ydim))
title('estimated spatial pattern')




 




