function result = ACfield_signal_lsq02(Bcoeffs, parms, data)
% function ACfield_signal_lsq02(Bcoeffs, parms, data)
%
% if called with two parameters only, this function generates an MR signal
% based on the signal equation with an extra term that accounts for the
% extra phase induced by an external oscillating magnetic field
%
% if a third input parameter is added (data), we assume that it's been
% called by a fitting routing, like lsqnonlin.  In that case, the function returns the RSS
% difference between those data and a signal generated from input parameters.
%
% all the input parameters are stuck into the structure 'parms'
%

gambar = 42.57e3;        % gamma/2pi in kHz/T
gam = gambar*2*pi;

%global parms
F = parms.F;    % FFT matrix
M0 = parms.M0;  % the object
ACfun = parms.ACfun;     % B field waveform (in time)
fovxy = parms.fovxy;

%Bcoeffs


% Generate space sampling location grid:

[xdim ydim] = size(M0);
kxdim = size(M0,1);
kydim = size(M0,2);

xmax = fovxy/2;
ymax = fovxy/2;
xstep = fovxy/size(M0,2);
ystep = fovxy/size(M0,1);
[x y] = meshgrid(-xmax:xstep: xmax-xstep  ,  -ymax:ystep: ymax-ystep );
x = x(:);
y = y(:);

% Build a magnetic field distribution from the DCT coefficients
ACfield = zeros(size(M0(:)));
Ncoeffs = length(Bcoeffs);
% for n=1:Ncoeffs
%     ACfield = ...
%         ACfield + ...
%         Bcoeffs(1,n)*cos((n-1)*x/xmax*pi) + ...
%         Bcoeffs(2,n)*sin((n-1)*x/xmax*pi) + ...
%         Bcoeffs(3,n)*cos((n-1)*y/ymax*pi) +...
%         Bcoeffs(4,n)*sin((n-1)*y/ymax*pi);
% end


% alternative code using polynomial bases
for n=1:Ncoeffs
    ACfield = ...
        ACfield + ...
        Bcoeffs(1,n)* (x/xmax).^(n-1) + ...
        Bcoeffs(2,n)* (y/ymax).^(n-1) ; 
end



% %ACfield = ACfield*1e-4;
 imagesc(reshape(ACfield,xdim,ydim)); 
 drawnow

ACfun2 = repmat(ACfun, kydim,1);
ACphase = gambar * cumsum( ACfun2,2 );
ACphase = ACphase(:)';

F_ac = exp(i * kron(ACfield(:), ACphase)); 

% Next, the dot-multiplication makes sure that each pixel's k-space trajectory includes the
% external field distorsion
% whos F F_ac F2
F2 = (F .* F_ac).';
signal2 = F2 * M0(:);


if nargin>2
    result = norm(double( abs (data -  signal2).^2)) ; % + (parms.beta * norm(gx(:)+gy(:)));
else
    result = signal2;
end


return

