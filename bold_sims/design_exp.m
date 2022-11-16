close all
points = 1;
iti = rand(50,1)*12 + 2;
iti = abs(iti);
iti(1) = iti(1) + 45;
iti(10) = iti(10) + 45;
iti(25) = iti(25) + 45;
iti(30) = iti(30) + 45;
times = cumsum(iti);

% create a set of stimuli occurring at random points in the time series
stim=zeros(720*points,1);
stim(round(times)) =1;




subplot 312, plot(stim),title('A set of stimuli at randomized times')

% create the HRF using SPM's function
hrf = spm_hrf(1/points);
subplot 311, plot(hrf),title('Response to a single stimulus')

% convolve the input and the HRF and clip out the end (empty)
resp=conv(stim,hrf);
resp = resp(1:size(stim));
subplot 313, plot(resp),title('Predicted Response to the stimuli')
figure
plot(abs(fftshift(fft(resp))))
