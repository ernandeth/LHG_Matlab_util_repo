clear
err = 5;
noise_level = 5;

NSAMPLES = 1200;
DURATION = 2400; % seconds
SECONDS = 100; % number of points in one second
NREGRESSORS = 5;
c = [ 1 0 0 0 0 ];
%iti = 4;
%noise_level = 3;
%tau = 1.5;

noise = noise_level *make1overf(NSAMPLES, 0.5);
fn = fft(noise);

beta = [1 2 3 4 5]';

% Create the Fake Design Matrix (model_resp):
model_resp = [];
model_error = [];

% create the HRFs

%hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
%hrf_error = make_hrf(2.5*SECONDS,(tau+err)*SECONDS, 20*SECONDS);

% Using the hrfs from SPM ...
p = [6 16 1 1 6 0 32 ];
hrf = spm_hrf( 1/SECONDS , p);
hrf = hrf/max(hrf);

p(1) = p(1) + err;
hrf_error = spm_hrf(1/SECONDS, p);
hrf_error = hrf_error/max(hrf_error);



plot(hrf)
hold
plot(hrf_error,'r')
hold
title('HRFs')
pause

p = [6 16 1 1 6 0 32 ];
h = spm_hrf( 2 , p);
h = h/max(h);

p(1) = p(1) + err;
he = spm_hrf( 2 , p);
he = he/max(he);

fh = abs(fftshift(fft(h)));
fhe = abs(fftshift(fft(he)));
fn = abs(fftshift(fft(noise)));

plot([-8:8],fh)
hold
plot([-8:8], fhe,'r');
plot([-599:600],0.01* fn, 'g');
%plot (fn(size(fn)/2 - size(fhe)/2 : size(fn)/2 + size(fhe)/2 ), 'g');
hold
title('fft of the HRFs and the noise')
%axis([0 1200 0 100])
pause

% Use a standard set of stimulus onsets instead ...
load 'fake_iti=4.mat'



% convolve the input and the HRF    
% The model response is our prediction (design matrix)
model_resp=conv2(stim,hrf);
model_resp = model_resp(1:DURATION*SECONDS, :);
%whos model_resp stim

model_error = conv2(stim, hrf_error);
model_error = model_error(1:DURATION*SECONDS, :);

% we will use the first regressor to explore what happens when 
% there is an unexpected time shift
shifting_reg = model_error(:,1);
% zero pad the regressor, so that it can be shifted in time
pad = zeros(1*SECONDS, 1);
shifting_reg = ...
   [pad;...
   shifting_reg;...
   pad];


% resample the response every TR
model_resp = model_resp(1: 2*SECONDS : DURATION*SECONDS,  : );
model_error = model_error(1: 2*SECONDS : DURATION*SECONDS,  : );

%whos model_resp model_error noise
%imagesc(model_resp)
%colormap(gray)
%pause

% Create the response data by weighting all the regressors by a beta parameter
% adding all the regressors together
% and adding noise to the result
data = model_resp * beta ; %+ noise;
ideal = model_error*beta ;


plot(fftshift(abs(fft(ideal))));
axis([1 1200 0 1000])
hold
pause

plot(fftshift(abs(fft(data))),'r');
pause

plot(fftshift(abs(fft(noise))), 'g');
title('fft of the data, the model and the noise')
hold
