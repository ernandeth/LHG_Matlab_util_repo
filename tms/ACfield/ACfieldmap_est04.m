function ACfieldmap_est04
% This is a simulation of the signal if we play a burst of B field between
% excitation and acquisition

%function ACfield_est = ACfieldmap_est(rawsignal, M0, Boff, fovxy, freq, delay, aqbw)

xdim=64;

M0 = phantom(xdim);
[xdim ydim] = size(M0);

freq = 1e3; % ACfield's frequency in Hz
rawsignal = fft2(M0);  % observed signal
fovxy = 20; % field of view
dt = 1e-6;

delay = 1e-4;

[xdim ydim] = size(M0);
[kxdim kydim] = size(rawsignal);

Slevel = max(abs(rawsignal(:)));

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
%Br = ACfield;

% A field that is NOT defined by a basis set:
Br = exp(-((x).^2 + (y-5).^2)/30);
Br = Br - exp(-((x).^2 + (y+5).^2)/30);


Br = Br / max(Br(:));

% test with uniform field:
% Br(:)=1e-8;

M0 = M0(:);
threshold = 0.02* mean(M0);
inds = find(abs(M0)>=threshold);


% start building  the parts of the signal equation integral:
% this is part with spin density times the the desired k-space trajectories:
% AND any DC offset in the field

F = zeros(length(M0(:)), length(rawsignal(:)));
Fi = zeros(length(M0(:)), length(rawsignal(:)));
for r=1:length(x)
    F(r,:) = exp(-j * (kx * x(r) + ky * y(r)) ) ;% .* exp (i * gambar * BoffDC(r) * t);
    Fi(r,:) = exp(j * (kx * x(r) + ky * y(r)) ) ;
end

% the field burst ends after 3 ms
maxt = 3.5e-3;
noiselevel = 0.01;
Bmax = 1e-5;
freq = 1e3; % ACfield's frequency in Hz


% These are the different parms that we will explore
Bmax_range = 10.^[-8.5:0.05:-3.5 ];
noise_range = 10.^[-6:0.1:1];
freq_range = 10.^[-1:0.1:5] ;
maxt_range = [8:0.05:12] * 1e-3;


%

count = 1;
NRMS1 = zeros(size(Bmax_range));

for Bmax = Bmax_range
    % the next part of the equation is the phase gained because of the AC field
    % this phase will be weighted by the field's amplitude at each pixel
    % when we update the AC field map
    tt = [0:dt:maxt];
    Br_scaled = Br * Bmax;
    
    ACfun = sin(2*pi*freq * tt);
    
    NRMS1(count) = computeField (tt, ACfun, Br_scaled, noiselevel, F, M0, xdim, dt);
    count = count+1;
end

figure (5)
subplot(221)
plot((Bmax_range), log(NRMS1))
xlabel('Max Field (Tesla)'); ylabel('Log NMRSE');

%

Bmax = 1e-5;
NRMS2 = zeros(size(noise_range));
count = 1;
for noiselevel = noise_range
    % the next part of the equation is the phase gained because of the AC field
    % this phase will be weighted by the field's amplitude at each pixel
    % when we update the AC field map
    tt = [0:dt:maxt];
    Br_scaled = Br * Bmax;
    
    ACfun = sin(2*pi*freq * tt);
    
    NRMS2(count) = computeField (tt, ACfun, Br_scaled, noiselevel, F, M0, xdim, dt);
    count = count+1;
end
figure (5)
subplot(222)
plot(log(Slevel./noise_range), log(NRMS2))
xlabel('SNR (dB)'); ylabel('Log NMRSE');

%

noiselevel = 0.05;

NRMS3 = zeros(size(freq_range));
count=1;
for freq = freq_range
    % the next part of the equation is the phase gained because of the AC field
    % this phase will be weighted by the field's amplitude at each pixel
    % when we update the AC field map

    tt = [0:dt:maxt];
    Br_scaled = Br * Bmax;
    
    ACfun = sin(2*pi*freq * tt);
    
    NRMS3(count) = computeField (tt, ACfun, Br_scaled, noiselevel, F, M0, xdim, dt);
    count = count+1;
end
figure (5)
subplot(223)
plot(freq_range, log(NRMS3))
xlabel('Signal Frequency (Hz)'); ylabel('Log NMRSE');

%

freq = 1e3; % ACfield's frequency in Hz

NRMS4 = zeros(size(maxt_range));
count=1;
for maxt = maxt_range
    % the next part of the equation is the phase gained because of the AC field
    % this phase will be weighted by the field's amplitude at each pixel
    % when we update the AC field map
    tt = [0:dt:maxt];
    Br_scaled = Br * Bmax;
    
    ACfun = sin(2*pi*freq * tt);
    
    NRMS4(count) = computeField (tt, ACfun, Br_scaled, noiselevel, F, M0, xdim, dt);
    count = count+1;
end
figure (5)
subplot(224)
plot(maxt_range, log(NRMS4))
xlabel('Burst Duration'); ylabel('Log NMRSE');

print -depsc ACfield_before_AQ
%

% THe next case simulates the Magstim coil:
tt=[0:dt:300e-6];
ACfun = -cos(3300*2*pi*tt) .* exp( -1.5e3 .* tt);
ACfun = [zeros(1,100) ACfun zeros(1,100)];
tt = [0:dt:500e-6];

Br_scaled = Br * Bmax;
NRMS = computeField (tt, ACfun, Br_scaled, noiselevel, F, M0, xdim, dt);

print -depsc Magstim_before_AQ

%
save ACfieldmap_est_wkspc

return
%%


function NRMS = computeField (tt, ACfun, Br, noiselevel, F, M0, xdim, dt)

gambar = 42.57e6;        % gamma/2pi in Hz/T
gam = gambar*2*pi;       % gamma in (rad/s) / T


doPlots = 1;

ACphase = gam * cumsum(ACfun)* dt ;
xtra_phase = exp(j*ACphase(end)* Br);

rawsignal = F * (M0(:).* xtra_phase(:));
rawsignal_clean =  F*M0(:);

noise = noiselevel *(randn(size(rawsignal)) + i*randn(size(rawsignal)));
rawsignal = rawsignal + noise;
noise = noiselevel *(randn(size(rawsignal)) + i*randn(size(rawsignal)));
rawsignal_clean = rawsignal_clean + noise;

im_distorted = fftshift(ifft2(fftshift(reshape(rawsignal,xdim,xdim))));
im_clean = fftshift(ifft2(fftshift(reshape(rawsignal_clean,xdim,xdim))));

phase_diff = -angle( im_distorted ./im_clean);

mask = zeros(size(M0));
mask(M0>0.1) = 1;
mask = reshape(mask,xdim,xdim);
phase_diff = phase_diff.*mask;

Br_est = phase_diff./ (eps+ACphase(end));
ACphase(end)

%e = (Br_est(:) - Br(:)) ./ Br(:) ;
e = (Br_est(:) - Br(:));
e = e.*mask(:);
NRMS = sqrt( e' *e )/ (xdim);
%NRMS = mean(abs(e)) ;

if doPlots
    figure(1),
    subplot(221)
    hold off; plot(tt, ACfun);
    %hold on; plot(tt, ACphase,'r');
    title('time function')
    
    subplot(222)
    imagesc(reshape(Br, xdim,xdim))
    title('spatial pattern')
    colorbar
    
    subplot(224)
    imagesc((abs(reshape(Br_est,xdim,xdim))))
    colorbar
    title('Estimated Field')
    
    subplot(223)
    imagesc((abs(reshape(log(e),xdim,xdim))))
    colorbar
    title('Log NMRSE')
    
    %{
    figure(2)
    imagesc(angle(F));
    title('FT matrix')
    figure(5)

    imagesc(phase_diff)
    colorbar; title('phase difference')
    drawnow
    
    %
    figure(4),
    subplot(211)
    imagesc(abs(reshape(rawsignal,xdim,xdim)));
    title('raw signal')
    
    figure(4), subplot(212)
    imagesc(abs(im_distorted ))
    title('distorted, reconned image')
    %}
    drawnow
end

return



