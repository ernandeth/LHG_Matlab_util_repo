% function ACfieldmap_gen()

M0 = phantom(32);

freq = 1e3; % ACfield's frequency
aqbw = 12e4;  % acquisition bandwith in Hz
delay = 0;

gambar = 42.57e3;        % gamma/2pi in kHz/T
gam = gambar*2*pi;

fovxy = 20; % field of view
xmax = fovxy/2;
ymax = fovxy/2;
xstep = fovxy/size(M0,2);
ystep = fovxy/size(M0,1);[x y] = meshgrid(-xmax:xstep: xmax-xstep  , -ymax:ystep: ymax-ystep );
[xdim ydim] = size(M0);

rawsignal = fft2(M0);  % observed signal
[kxdim kydim] = size(rawsignal);
rawsignal = rawsignal(:)';

% define the spatial shape of the DC and the AC fields
BoffDC = 0 * exp(-(x.^2 + y.^2)/5);
BoffAC = 10 * exp(-((x-2).^2 + (y-2).^2)/5);


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


% shortcut:  figure out which pixels have enough signal to contribute
M0 = M0(:);
threshold = 0.02* mean(M0);
inds = find(abs(M0)>=threshold);

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
% the resulting matrix has dimensions of N_pixels by N_timesamples
t = linspace(0, kxdim/aqbw, kxdim);
t = repmat(t, kydim,1);
ACphase = gambar * cumsum( sin(2*pi/aqbw*(freq*t + delay)),2 );
ACphase = ACphase(:)';
ACphase_r = kron(BoffAC(:), ACphase);

% Note, the kron simply makes sure that each pixel's contribution to the
% net phase of the signal is weighted appropriately

% now include that AC field phase into the contribution of that pixel to the signal
% and compute the integral
signal = sum(M0phase .* exp( i * ACphase_r ) ,1 );
signal =reshape(signal,xdim,xdim);
im1 = fftshift(fft2(fftshift(signal)));


%% method 2: the linear algebra method:
% construct the Fourier Matrix
F = zeros(length(M0), length(rawsignal), 'single');
for r=1:length(x)
    F(r,:) = exp(i * (kx * x(r) + ky * y(r)) ) ;% .* exp (i * gambar * BoffDC(r) * t);
end


% construct the distortion matrix produced by AC field (F_ac):
% Note, the kron product makes sure that each pixel's contribution to the
% net phase of the signal is weighted appropriately
F_ac = exp(i * kron(BoffAC(:), ACphase)); 

% Next, the dot-multiplication makes sure that each pixel's k-space trajectory includes the
% external field distorsion
F2 = (F .* F_ac)';

figure(1), subplot(211)
imagesc(angle(F2));

signal2 = F2 * M0;
signal2 =reshape(signal2,xdim,xdim);

% something went weird and we need to flip the data:
signal2 = signal2(end:-1:1,end:-1:1);

im2 = fftshift(fft2(fftshift(signal2)));

%% Solving for BoffAC : collect second image with phase shift
% Calculate a second signal with a phase shift in the AC waveform:
% caculate a signal in the presence of the oscillating field (i.e., the
% integral of the waveform)
delay = 1e-2;
t = linspace(0, kxdim/aqbw, kxdim);
t = repmat(t, kydim,1);
ACphase_b = gambar * cumsum( sin(2*pi/aqbw*(freq*t + delay)),2 );
ACphase_b = ACphase_b(:)';

ACphase_r = kron(BoffAC(:), ACphase_b);

% construct the distortion matrix produced by AC field (F_ac):
% Note, the kron product makes sure that each pixel's contribution to the
% net phase of the signal is weighted appropriately
F_ac = exp(i * kron(BoffAC(:), ACphase_b)); 

% Next, the dot-multiplication makes sure that each pixel's k-space trajectory includes the
% external field distorsion
F2b = (F .* F_ac)';

figure(1), subplot(212)
imagesc(angle(F2b));

signal2b = F2b * M0;
signal2b =reshape(signal2b,xdim,xdim);

% something went weird and we need to flip the data:
signal2b = signal2b(end:-1:1,end:-1:1);

im2b = fftshift(fft2(fftshift(signal2b)));

delta_phs = angle(im2 ./ im2b); 
Bsolve = delta_phs / (1 - exp(i*2*pi*delay/aqbw));
Bsolve2 = fftshift(fft2(fftshift(Bsolve)));

%% compare the signals and the distorted images.
figure(2)
subplot(231)
imagesc(angle(signal))
subplot(232)
imagesc(angle(signal2))
subplot(233)
imagesc(angle(signal2b))

subplot(234)
imagesc(abs(im1));
subplot(235)
imagesc(angle(im2));
subplot(236)
imagesc(angle(im2b));

figure(3)
subplot(211)
imagesc(delta_phs);
subplot(212)
imagesc(abs(BoffAC))





% return
