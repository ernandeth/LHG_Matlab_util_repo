close all

cutoff = 0.15;
TR = .32;
nyquist = 1/(2*TR);
N_coeffs = 10;

%%%  Scott : change this:
%data = sin(0:0.01:2*pi);
data = card/max(data);
fdata =fftshift(fft(data));
%%%%%%%%%%%%%%%%%%%%%%


f = linspace(-nyquist, nyquist, size(data,2)+1);
f = f(1:end-1);

% convert cutoff frequency from Hz to pi*rads
cutoff = cutoff  / nyquist;
b = remez(N_coeffs, [0 cutoff-0.05  cutoff+0.05  1], [1 1 0 0]);

impulse= zeros(size(data));
impulse(max(size(data))/2) = 1;

impulse_response = filtfilt(b,1,impulse);
freq_response = fftshift(fft(impulse_response));
freq_impulse = fftshift(fft(impulse));

freq_gain = 20*log10( abs(freq_response).^2 ./ abs(freq_impulse) ) ;

filtered = filtfilt(b,1,data);
filtered = filtered/max(filtered);

subplot(321), hold on, plot(f, freq_gain , 'g')
plot(f,abs(fftshift(fft(data)))-200,'k');
title('using remez filter');
axis tight;

subplot(322), plot(data,'k');
hold on
subplot(322), plot( filtered,'g');
axis tight;

% 
rem_filt = freq_response;

%%%  Now plot the frequency resplonse of multi-shot imaging
z = exp(-j*pi*f/nyquist);

% Two shots:
%freq_response = (sqrt(2 + 2*cos(f*pi/nyquist)));
freq_response = (1 + z).* z.^-0.5 ;
gain_dB = 20*log10(abs(freq_response));


subplot(323), plot(f,gain_dB,'b')
hold on,
plot(f,100*angle(freq_response))
title('simple averaging 2 shots')
axis tight

% apply the filter
tmp = freq_response .* fdata .*rem_filt;
filtered = abs(fftshift(fft(tmp)));
filtered = filtered/max(filtered);

subplot(324), plot(data,'k');
hold on
subplot(324), plot(filtered,'b');
axis tight

filt_ratio2 = abs(rem_filt) ./abs(freq_response);

% four shots:

freq_response = (1 + z +  z.^2 + z.^3).* z.^-1.5 ;
gain_dB = 20*log10(abs(freq_response));

subplot(325), plot(f,gain_dB,'r')
axis([-nyquist nyquist -100 10])
title('simple averaging four shots')

% apply the filter
tmp = freq_response .* fdata;
filtered = abs(fftshift(fft(tmp)));
filtered = filtered/max(filtered);

subplot(326), plot(data,'k'), hold on;
hold on
subplot(326), plot(filtered,'r');
axis tight


filt_ratio = abs(rem_filt)./abs(freq_response);
figure
plot(abs(filt_ratio));
figure
plot(abs(filt_ratio2));

