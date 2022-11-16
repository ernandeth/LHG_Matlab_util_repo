function t_score = resp_gen_ana(noise_level, iti)
% function t_score = resp_gen_ana(noise_level,  iti)
%
% Luis Hernandez
% University of Michigan
%
% This program predicts the t_score of the 
% correlation coefficient between
%  - a simulated BOLD response data set, and ..
%  - a model of those data that is time shifted relative to the data
% as a function of that time shift, and the noise.  
% The correlation coefficient is computed 
% using an ANALYTICAL expression, rather than simulating data and 
% computing the correlation coefficients for each time shift.
%
% Arguments:
% noise_level - Amplitude of the smallest amount of noise to be added.
%          this will be incremented up to 10 times to show a family of curves.
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
%- The canonical HRF used is the one that SPM uses spm_hrf()
%- The experiment is assumed to last DURATION seconds
%- TR = 2seconds
%

SECONDS = 100; % number of points in one second
DURATION = 600;

% create a set of stimuli occurring at random points in the time series
stim=zeros(DURATION*SECONDS,1);
times = (1:iti:DURATION);
%times = SECONDS*(times + (iti/2)*rand

% Getting a Gaussian Distribution of ITI's
mean_iti = iti;
iti = (mean_iti/2)*randn(size(times));
iti = iti + abs(min(iti));
times = times + iti;

times = times*SECONDS;
stim(times) =1;


% create the HRF
%hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
hrf = spm_hrf(1/SECONDS);
hrf = hrf/max(hrf);


disp('convolving the input and the HRF and clip out the end')
resp=conv(stim,hrf);
resp = resp(1:SECONDS:DURATION*SECONDS);
disp('done')

plot(resp), title('Response')
pause
% Compute the correlation coefficient between the delayed data and the model
% The derived analytical expression is
%             { auto_corr(delay) - mean(resp)^2 }
% rho =       --------------------------------------
%         {std(resp)*sqrt( std(resp)^2 + std(noise)^2) }


% Compute the autocorrelation from -1 second every 100 ms.until 1 second
auto_corr = xcorr(resp,'unbiased');
midpoint =size(auto_corr)/2;
ac = auto_corr(midpoint - 1*SECONDS: 0.1*SECONDS: midpoint + 1*SECONDS)'
mean_resp = mean(resp)
std_resp = sqrt(var(resp))

% auto_corr = auto_corr(size(auto_corr)/2 : size(auto_corr));
% subplot 311, plot(auto_corr,'r'), title ('autocorrelation');
% subplot 312, plot(abs(fft(auto_corr))), title('Power Spectrum');

% increment the noise level from 0 to 30 times the noise_level
for i=1:10
   std_noise = i * noise_level;
   rho = [rho;
      (ac - mean_resp^2 )./(std_resp * sqrt(std_resp^2 + std_noise^2)) ];
end
t_score = rho;
%t_score=rho * sqrt(1198) ./ sqrt(1 - rho.^2);

%subplot 312, imagesc(rho)
%colorbar

return