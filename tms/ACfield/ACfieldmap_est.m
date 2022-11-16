%function ACfield_est = ACfieldmap_est(rawsignal, M0, Boff, fovxy, freq, delay, aqbw)

DEBUG=1;
if DEBUG
	M0 = phantom(32);
	
	freq = 1e3; % ACfield's frequency
	aqbw = 12e3;  % acquisition bandwith in Hz
	rawsignal = fft2(M0);  % observed signal
	fovxy = 20; % field of view


	delay = 1e-4;
	xmax = fovxy/2;
	ymax = fovxy/2;
	xstep = fovxy/size(M0,2);
	ystep = fovxy/size(M0,1);[x y] = meshgrid(-xmax:xstep: xmax-xstep  , -ymax:ystep: ymax-ystep );

	BoffDC = 0.01 * exp(-(x.^2 + y.^2)/5);
	BoffAC = 0.01*exp(-((x-2).^2 + (y-2).^2)/5);
	
end

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

% if isSpiral
%	% kx, ky = genspiral ....
%	kdata = griddata(  k(:,1), k(:,2), rawsignal,kxc, kyc);
%	kdata(find(isnan(kdata))) = 0;
% end

M0 = M0(:);
threshold = 0.02* mean(M0);
inds = find(abs(M0)>=threshold);

rawsignal = rawsignal(:)';

t = linspace(0, kxdim/aqbw, kxdim);
t = repmat(t, 1, kydim);

M0phase = zeros(length(M0), length(rawsignal), 'single');
% start building  the parts of the signal equation integral:
% this is part with spin density times the the desired k-space trajectories:
% AND any DC offset in the field

for r=1:length(x)
	fprintf('\rContribution from pixel  %d ', r);
	M0phase(r,:) = M0(r) * ...
		exp(i * (kx * x(r) + ky * y(r)) ) .* ...
		exp (i * gambar * BoffDC(r) * t);
end

% the next part of the equation is the phase gained because of the AC field
% this phase will be weighted by the field's amplitude at each pixel
% when we update the AC field map
ACphase = gambar * cumsum( sin(2*pi*(freq/aqbw)*t + delay/aqbw) );
ACphase = single(ACphase);
ACphase = repmat(ACphase, 1, size(rawsignal,1));

parms.M0phase = M0phase;
parms.ACphase = ACphase;
parms.inds = inds;

dummy = [];

if DEBUG
    tstart = tic;
    
	est = ACfield_signal_lsq(BoffAC, parms); 
    
    fprintf('\nFirst call to the Integrating function takes %f secs\n',toc(tstart));
    
	est = reshape(est,xdim,ydim);
	subplot(211)
	imagesc(abs(est))
	f = fftshift(fft2(fftshift(est)));
	subplot(212)
	imagesc(abs(f));
    
    rawsignal = est;
	drawnow
end

optvar=optimset('lsqnonlin');
optvar.Display = 'iter';
optvar.MaxIter = 50;
optvar.MaxFunEvals = (length(rawsignal(:))+1)*400;
optvar.TolFun = 1e-10;

% initial estimate of the AC field: there is none 
est0  = zeros(size(M0));
%est0 = BoffDC;



tstart = tic;
% call the fitting routine:
[ACfield_est res] = lsqnonlin(@ACfield_signal_lsq, est0,...
	[],[],...
	optvar,...
	parms,...
	rawsignal(:)');

fprintf('\nTime Elapsed = %f',toc(tstart));
figure
imagesc(reshape(ACfield_est,xdim,xdim))

% return
