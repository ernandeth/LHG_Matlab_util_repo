function t_score = bold_shift_ex(noise_level, tau, iti)
% function t_score = bold_shift_ex(noise_level, tau, iti)
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
% noise_level - std deviation of the  noise to be added.
%          the noise is expressed as the fraction of the BOLD amplitude.
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

DURATION = 2400; % seconds
SECONDS = 10; % number of points in one second

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
resp=conv(stim,hrf);
resp = resp(1:DURATION*SECONDS);
%subplot 313, plot(resp), title('Response')
%plot(resp), title('Response')

% resample the response every 2 sec (TR)
% model_response = resp(1:2*SECONDS:size(resp));


% Compute the correlation coefficient between the delayed data and the model
% The derived analytical expression is
%             { auto_correlation(resp) - mean(resp)^2 }
% rho =       --------------------------------------
%         {std(resp)*sqrt( std(resp)^2 + std(noise)^2) }


auto_corr = xcorr(resp,randn(2400*SECONDS,1),'coeff');
xx = [-2399:2400];
auto_corr = auto_corr(2 : SECONDS: 4800*SECONDS -1);
whos
plot(xx, auto_corr)

%auto_corr = auto_corr(size(auto_corr)/2 : size(auto_corr));
%subplot 211, plot(auto_corr,'r'), title ('autocorrelation');
%subplot 212, plot(abs(fft(auto_corr))), title('Power Spectrum');

% resample the autocorrelation every 100 ms. from -1 sec. to 1 sec.
ac = xcorr(resp,1 * SECONDS, 'unbiased');
ac = ac(1: SECONDS/10: size(ac));

mean_resp = mean(resp);
std_resp = sqrt(var(resp));

% apply the equation here:
correlations = (ac - mean_resp^2 )./(std_resp * sqrt(std_resp^2 + noise_level)) ;

% convert correlation coefficent to t score:
df = size(resp,1)/ (SECONDS*2) -2;
whos
t_score = correlations * sqrt(df) ./ sqrt(1 - correlations.^2)

figure
shift = [-1:0.1:1];
plot(shift, t_score)



return

