function correlations = resp_gen_rand(noise_level, tau, iti)
% function correlations = resp_gen_rand(noise_level)
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM, 
% 
% Gaussian noise is added determined by 'noise_level' 
% (fraction of a single HRF)
%
% The function is then shifted in time for five different
% delays
%
% The total imaging time is 2400 sec.
%
% The number of points in one second is specified by the 
% variable 'points'.( not related to TR )
%

points = 100; % number of points in one second

% create a set of 50 stimuli occurring at random points in the time series
% the stimuli must be between 8 and 16 seconds apart.
stim=zeros(2400*points,1);
times = (1:iti:2400);
times = points*(times + 4*rand(size(times)));
stim(times) =1;
%subplot 311, plot(stim);

% create the HRF
hrf = make_hrf(2.5*points,tau*points, 20*points);
%hrf = make_hrf(2.5*points,0.1*points, 20*points);
%subplot 312, plot(hrf)

% convolve the input and the HRF and clip out the end (empty)
resp=conv(stim,hrf);
resp = resp(1:2400*points);

% The model response is what we predict
model_response = resp;

% add gaussian noise to the response data, and resample every TR.  
% This is the 'sampled data'
resp = resp + noise_level *randn(size(resp));
resp = resp(1:2*points:size(resp));

% shift the model by a series of delays: 100-2000 msec.
% and compute the correlation coefficient between the data and 
% the model


for delay=0: 0.1 : 1.9
   % resample the model every TR and shift it
   resp_d = model_response(1+ delay*points: 2*points: size(model_response));
   
   delay
   whos 
   
   p = corrcoef(resp, resp_d );
   correlations = [correlations   p(1,2)];
end
% we cheat on the last one ...
last = size(correlations,2);
correlations = [correlations correlations(:,last)];

% save sim_noise.mat resp_delays model_response correlations

return

%%%%%
function hrf = make_hrf(delta,tau, npoints)
% Generates a gamma variate function for the hemodynamic response function.
% points is the number of points to be included in the function.
% (ie - the length of the response to be calculated)
% Normalized so that it peaks at 1

   hrf=zeros(npoints,1);

	for t=1:size(hrf,1)
   		hrf(t)=((t-delta)./tau).^2.*exp(-(t-delta)./tau).*((t-delta) > 0);
   end  
   
   hrf = hrf/max(hrf);

return

%%%%
function result = make_squares(onsets, lengths)
	t_onset = onsets;
	t_duration = lengths;

	result=zeros(200 ,1);
	for i=1:size(t_onset,1)
   		result(floor(t_onset(i)):floor(t_onset(i) + t_duration(i)))=1;
	end
   
return
