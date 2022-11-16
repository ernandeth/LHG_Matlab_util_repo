function correlations = resp_gen_rand(noise_level)
% function correlations = resp_gen_rand(noise_level)
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM, 
% 
% Gaussian noise is added determined by 'noise_level' 
% (% of a single HRF)
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

% create a set of stimuli occurring at random points in the time series
stim=zeros(100*points,1);
stim(round(size(stim,1)*rand(10,1))) =1;

subplot 312, plot(stim),title('A set of stimuli at randomized times')
subplot 312, axis ([0 size(stim,1) -0.1 1.2])
axis off
pause

% creat the HRF
%hrf = make_hrf(2.5*points,1.25*points, 20*points);
hrf = spm_hrf(1/points);
subplot 311, plot(hrf),title('Response to a single stimulus')
axis off
pause

% convolve the input and the HRF and clip out the end (empty)
resp=conv(stim,hrf);
resp = resp(1:size(stim));
subplot 313, plot(resp),title('Predicted Response to the stimuli')
axis off
pause
return


% resample the response, so that we can use it for the 
% model to fit.
model_response = resp(1:2*points:size(resp));
subplot 313, plot(model_response),title('Predicted observation of the response ')
%axis([0 50 -1.5 1.5])
axis off
pause

% add gaussian noise to the response data
model_response = model_response + noise_level *randn(size(model_response));

hold on
subplot 313, plot(model_response,'r')

hold off

% shift the data by a series of delays: 100 - 500 msec.
%for delay=0: 0.2 : 1
%  resample the data to the acquisition rate and shift it
%	stim2 = stim(1+ delay*points: 2*points: size(stim));
%	hrf2 = hrf  (1+ delay*points: 2*points: size(hrf));
%	resp2 = resp(1+ delay*points: 2*points: size(resp));
%   
%   resp_delays = [resp_delays resp2];
%   
%end

%subplot 311, plot(hrf2)
%subplot 312, plot(stim)
%subplot 313, plot(resp2)

% compute the correlation coefficient between the delayed data and 
% the model

%for i=1:size(resp_delays,2);
%   p = corrcoef(model_response, resp_delays(:,i));
%   correlations = [correlations   p(1,2)];
%end


%save sim_noise.mat resp_delays model_response correlations

return

%%%%%
function hrf = make_hrf(delta,tau, npoints)
	   
   hrf=zeros(npoints,1);

	for t=1:size(hrf,1)
   		hrf(t)=((t-delta)./tau).^2.*exp(-(t-delta)./tau).*((t-delta) > 0);
   end  
      
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
