
noise_level = 0
tau = 6
iti = 4

DURATION = 2000 % seconds
points = 100; % number of points in one second

% create a set of stimuli occurring at random points in the time series
% the stimuli must be between 8 and 16 seconds apart.
stim=zeros(DURATION*points,1);
%times = (1:iti:DURATION);
%times = points*(times + 4*rand(size(times)));




% Actual regressors: (makes the iti parameter irrelevant)
load stim_onsets
times = stim_onsets * 2 * points;


stim(times) =1;
subplot 311, plot(stim);


% create the HRF
%hrf = make_hrf(2.5*points,tau*points, 20*points);
hrf = spm_hrf(0.01);

subplot 312, plot(hrf)


% convolve the input and the HRF and clip out the end (empty)
resp=conv(stim,hrf);
resp = resp(1:DURATION*points);
subplot 313, plot(resp), title('Response')
%plot(resp), title('Response')
pause
figure

% resample the response every 2 sec (TR)
% model_response = resp(1:2*points:size(resp));


% Compute the correlation coefficient between the delayed data and the model

% The derived analytical expression is
%             { auto_correlation(resp) - mean(resp)^2 }
% rho =       --------------------------------------
%         {std(resp)*sqrt( std(resp)^2 + std(noise)^2) }


%auto_corr = xcorr(resp,'unbiased');
%auto_corr = auto_corr(size(auto_corr)/2 : size(auto_corr));
%subplot 211, plot(auto_corr,'r'), title ('autocorrelation');
%subplot 212, plot(abs(fft(auto_corr))), title('Power Spectrum');

% resample the autocorrelation every 100 ms. from -1 sec. to 1 sec.
%ac = auto_corr(size(auto_corr)/2 - 1*points +1 : 0.1*points: size(auto_corr)/2 + 1*points +1)'
ac = xcorr(resp,0.5 * points, 'unbiased');
ac = ac(1:10:size(ac));

mean_resp = mean(resp)
std_resp = sqrt(var(resp))

% apply the equation here:
correlations = (ac - mean_resp^2 )./(std_resp * sqrt(std_resp^2 + noise_level^2)) ;

% convert correlation coefficent to t score:
t_score = correlations * sqrt(1198) ./ sqrt(1 - correlations.^2)

shift = [-0.5*points+1:100:0.5*points];
plot(shift, ac)
pause
plot(shift, t_score)

