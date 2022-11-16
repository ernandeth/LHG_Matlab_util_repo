function t_score = bold_shift_ana10(noise_level, tau, iti)
% function t_score = bold_shift_ana10(noise_level, tau, iti)
%
% Luis Hernandez
% University of Michigan
%
% This program predicts the correlation coefficient between
%  - a simulated BOLD response data set, and ..
%  - a model of those data that is time shifted relative to the data
% as a function of that time shift, and the noise.  
% The correlation coefficient is computed 
% using an ANALYTICAL expression, rather than simulating data and 
% computing the correlation coefficients for each time shift.  
%
% The function returns the corresponding T-score of the correlation coefficient
%
% Arguments:
% noise_level - variance of the smallest amount of noise to be added.
%          the variance of the noise is expressed as 
%          a fraction of the BOLD amplitude.
%          the variance is stepped up to 10 times to generate
%          a familiy of curves.
% tau  -   this parameter determines the width of the canonical HRF
% iti   -  mean inter-trial-interval
%
% the function generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM, 
% 
% other notes:
% - The total imaging time is DURATION sec.
%
% - The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%

DURATION = 600; % seconds
SECONDS = 10; % number of points in one second
TR =1.0;

% create a set of stimuli occurring at random points in the time series
% the stimuli must be between 8 and 16 seconds apart.
stim=zeros(DURATION*SECONDS,1);
times = (1:iti:DURATION);
times = SECONDS*(times + (iti/2)*rand(size(times)));
stim(times) =1;

% Actual regressors: (makes the iti parameter irrelevant)
%load stim_onsets
%times = stim_onsets * 2 * SECONDS;
%stim(times) =1;

%subplot 311, plot(stim);


% create the HRF using SPM's function
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

%subplot 312, plot(hrf)


% convolve the input and the HRF and clip out the end (empty)
resp = conv(stim,hrf);
resp = resp(1: DURATION*SECONDS);

%subplot 313, plot(resp), title('Response')
%

% resample the response every 2 sec (TR)
% model_response = resp(1:2*SECONDS:size(resp));


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
ac = xcorr(resp, 1*SECONDS, 'unbiased'); 
%figure
%subplot 211, plot(ac);

ac = ac(1:0.1*SECONDS:end);

mean_resp = mean(resp)
std_resp = sqrt(var(resp))

% apply the equation here:
rho = [];
for i = 1:10
   
   noise = i*noise_level;
   
   correlations = (ac - mean_resp^2 ) ./ (std_resp * sqrt(std_resp^2 + noise)) ;
   
   rho = [rho  correlations];
end

% convert correlation coefficent to t score:
df = size(resp,1)/ (TR*SECONDS) - 2
t_score = rho .* sqrt(df) ./ sqrt(1 - rho.^2)

%t_score = rho;

%subplot 212, plot (t_score)

%shift = [-0.5*SECONDS+1:100:0.5*SECONDS];

return

