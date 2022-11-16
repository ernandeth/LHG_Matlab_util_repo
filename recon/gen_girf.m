function [realp imagp maxFreq] =  gen_girf(fid_file)
% function [realp_coeffs imagp_coeffs maxFreq_inChirp] =  gen_girf(fid_file)
%
% calculate gradient input response function from a chirp measurement
% fid_file is the .fid directory for the CHirp data (contains .fid and
% propcpar)
%
% e.g.:
% chirp_05.fid: coronal -> z
% chirp_06.fid : axial  -> y
% chirp_07.fid : sagittal -> x
%
% returns the polynomial coefficients of the real and imaginary components
% of the transfer function as well as the max. Frequency of the transfer
% function in Hz.
% 
gamma = 4257.748;  % Hz / Gauss

Gis = [];
nomGrad = [];
realp = [];
imagp = [];
H_nice = [];
H2_nice = [];


% Import and process data
[ A procpar] = getFid(fid_file);

% only the first 10 K points are good... too noisy afterward
%     NPOINTS = 5e3;
%     A = A(1:NPOINTS,:);

[a b] = size(A);

% every other echo has a gradient trajectory turned on
% first half of signals = reference FID
% second half = FID in presence of gradient
idx1 = 1:b/2;
idx2 = b/2+1:b;

dt = 1/procpar.sw;
pos = procpar.pss;
phiref = A(:,idx2);
phiRaw = A(:,idx1);

%%%%% error in z axis: there are only four positions, not six %%%%%
pos = pos(1:4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subtract reference phase
tmp = (phiRaw ./ phiref );
phi = zeros(size(phiref));

for j=1:4
    phi(:,j) = phase(tmp(:,j));
end


Gis = zeros(length(phi) , length(pos));
count = 1;
% for every slice collected we have a separate measurement of the
% trajectory
for m = 1:length(pos)
    
    % calculate the derivative in time at each time point
    for p = 2:length(phi)-1
        
        Gis(p, count) = ( phi(p+1, m) - phi(p-1, m)) / (gamma*dt*4*pi*pos(m));
    end
    
    Gis(1, count) = (phi(2,m) - phi(1,m)) / (gamma*dt*4*pi*pos(m));
    Gis(end, count) = (phi(end,m) - phi(end-1,m)) / (gamma*dt*4*pi*pos(m));
    count = count + 1;
end
allpos = Gis;

NPOINTS = length(Gis);

% (average all the slices) we will have ONE gradient waveform per axis
Gis = mean(Gis,2);

% Get th nominal gradient wavefor:  a 'chirp' function :
[inputChirp tsamp] = getInputChirp(procpar);
inputChirp = inputChirp(1:NPOINTS);

nomGrad = [nomGrad  inputChirp'];

% There is a 40 ms. delay in the gradient waveform
% throw out the first 40 ms.
% hence, we don't get the last 40 ms of the nominal Chirp either
% aqdelay = 40e-6 / dt;
aqdelay = 0;  % or ... may be not.

Gis = Gis(aqdelay+1:end , :);
nomGrad = nomGrad(1 : end-aqdelay, :);

% Now get rid of another 100 points because of ring-down from prephasing
% gradient.
Gis = Gis(4*aqdelay+1 : end , :);
nomGrad = nomGrad(4*aqdelay+1 : end, :);

% calculate frequency domain of input
nomGradF = fft(nomGrad,[],1);

% calculate frequency content of measured output
GisF = fft(Gis,[],1);

% Calculate transfer function in frequency domain
H = GisF ./ (eps + nomGradF);


% a time vector
Npoints = size(Gis,1);

t = [0:Npoints-1]*dt;
% Nyquist frequency for the input waveform
% Note that this is not the maximum frequency in the chirp/
Nyq = procpar.sw/2;
df = 2*Nyq/Npoints;
f = linspace(-Nyq, Nyq-df, Npoints);



% we;re only interested in the stuff we played.
% cutoff frequency is the max frequency of the input chirp
cutoff = procpar.f1;

cutoff = (Nyq-cutoff)/2;
cutoff = Nyq;

cutoff_N = floor(cutoff * (Npoints/2 -1) / Nyq);
maxFreq = cutoff;

% f2 and H2 are the lower portion of the spectrum:
% where our input chirp is
f2 = linspace(-cutoff, cutoff-df, cutoff_N*2);
f2pos = f2(end/2:end);

H2 = [...
    H(1:cutoff_N,:) ;
    H(end-cutoff_N+1 : end, :) ];



% Show the original data in time and freq. domains
figure (1)
subplot(211)
plot(t, Gis)
hold on
plot(t, nomGrad,'r')
hold off
axis tight
title('Gradient waveforms ')
legend( 'Measured Output', 'Nominal Input')
hold off

subplot(212)
plot(f, fftshift(real(GisF)))
hold on
plot(f, fftshift(real(nomGradF)),'r')
hold off
axis tight
title('Gradient FFT (real) ')
legend('Measured Output' , 'Nominal Input')

hold off

drawnow


% Now try to fit a polynomial to each transfer function:
% real and imaginary separately
% Real part:
tmp = real(H2)  ;

% take the positive part of the spectrum only
tmp = tmp(1:end/2);

[p,s] = polyfit(linspace(0, cutoff-df ,length(tmp)), tmp', 50);
realjunk = polyval(p,linspace(0, cutoff-df ,length(tmp)));
realp = p;

% include a rampdown beyond the cutoff frequency
% rdown = tmp(end) * hanning(length(tmp)/4);
% tmp = [tmp; rdown(end/2+1:end)];
% [p,s] = polyfit(linspace(0, 9/8*cutoff-df ,length(tmp)), tmp', 20);


% rebuild the negative part of the spectrum
% remember:  I havent' done any fftshifts
realjunk(1) = 1;
realjunk = [ realjunk  0 realjunk(end:-1:2)];


%  Imaginary part:
tmp = imag(H2)  ;
tmp = tmp(1:end/2);


[p,s] = polyfit(linspace(0, cutoff-df ,length(tmp)), tmp', 50);
imjunk = polyval(p,linspace(0,cutoff-df, length(tmp)));
imagp = p;

% rebuild the negative part of the spectrum
% remember:  I havent' done any fftshifts
imjunk(1) =0;
imjunk = [imjunk 0 -imjunk(end:-1:2) ];


% the 'ideal' transfer function:
H2_nice = complex(realjunk, imjunk);

tmp = H2  ;

figure (2)
subplot(211)
plot(f2, fftshift(real(tmp)),'k')
hold on
plot(f2, fftshift(real(H2_nice)))
xlabel('Frequency (Hz)')
hold off
%axis([0 2000 0 1])
title('Real of Transfer Function')
legend('Measured', 'Estimated')

subplot(212)
plot(f2, fftshift(imag(tmp)),'k');
hold on
plot(f2, fftshift(imag(H2_nice)))
xlabel('Frequency (Hz)')
hold off
title('Imag of Transfer Function')
legend('Measured', 'Estimated')

figure (8)
subplot(211)
plot(f2, fftshift(abs(tmp)),'k')
hold on
plot(f2, fftshift(abs(H2_nice)))
xlabel('Frequency (Hz)')
hold off
%axis([0 2000 0 1])
title('mag of Transfer Function')
legend('Measured', 'Estimated')

subplot(212)
plot(f2, fftshift(angle(tmp)),'k');
hold on
plot(f2, fftshift(angle(H2_nice)))
xlabel('Frequency (Hz)')
hold off
title('Phase of Transfer Function')
legend('Measured', 'Estimated')


% upsample the Transfer function to the original resolution:
zpad = zeros(1, Npoints-length(H2_nice));
H_nice = [H2_nice(1:end/2) zpad H2_nice(end/2+1 : end) ];

% the gradient's impulse response function in the imte dpmain
% is the iFFT of the transfer function in frequency
girf = ifft(H_nice);
girf = girf(1:Npoints/2-1);


% try to predict the output given the input and the impulse response
% function:
output_nice = conv(nomGrad.', girf) ;
output_nice = output_nice(1:Npoints) ;

% Plot outputs and inputs.
t = linspace(0, procpar.at, Npoints)';
tmp = Gis;

figure (3)
subplot(211)
plot(t, real(output_nice),'k')
hold on
plot(t, real(tmp))
hold off
axis([0 0.05 -1 1])
title('Real of output  ' )
legend('Predicted', 'Measured ')

subplot(212)
plot(t, imag(output_nice),'k')
hold on
plot(t, imag(tmp))
title('Imaginary of output')
legend('Predicted', 'Measured ')
axis([0 0.05 -1e-3 1e-3])
hold off
drawnow


return
% save /Users/hernan/data/girf/girf_parms realp imagp H_nice dt Nyq Gisirf cutoff



function ft_girf(H)
% Doug's version using a Hanning window and Fourier basis
Hsh = fftshift(H);

subplot(211)
plot(abs(Hsh))
subplot(212)
plot(angle(Hsh))

csz = 1432; % part of spectrum to fit
cind = [-csz/2:csz/2-1]';
bind = [-3500:3499]';

% keep the interesting part
Hcr = Hsh(3500-csz/2 : 3500+csz/2-1);
% phase gain per sample
linpersamp = angle(sum(Hcr(1:end-1) .* conj(Hcr(2:end))));
% remove the phase gain
Hcrnp = Hcr.*exp(1i*linpersamp*cind);

subplot(211)
hold on, plot(abs(Hcrnp),'r'), hold off
subplot(212)
hold on, plot(angle(Hcrnp),'r'), hold off

% make hanning window for roll off at high frequency
f1 = hanning(csz*2);
zfill = zeros([(7000-3*csz)/2 1]);

filt = [zfill; 
    f1(1:end/2); 
    ones([csz 1]); 
    f1(end/2+1:end); 
    zfill];

plot(filt)

% apply filter and put the linear phase  back in.

Hbig = [zfill; 
    Hcrnp; 
    Hcrnp; 
    Hcrnp; 
    zfill].*filt.*exp(-1i*linpersamp*bind);

subplot(211)
plot(abs(Hbig))
subplot(212)
plot(angle(Hbig))

% Use these many fourier coefficients
nvals = 30;

% make another filter for the reduced transfer function
f1 = hanning(7000-csz)';
filt = [f1(1:end/2) ones([1 csz]) f1(end/2+1:end)];

mft = ifft(((Hbig)));

mft2 = zeros([1 7000]);
mft2(1:nvals) = mft(1:nvals);

mft2(end-nvals+1:end) = mft(end-nvals+1:end);
mftfit = abs(fft(mft2));
return
