function t_score = bold_shift_sim1(noise, tau, iti)
% function t_score = bold_shift_simi1(noise, tau, iti)
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM ITIs (with a Gaussian distribution), 
% It then computes the t-score of a General Linear Model REgression between the 
% model and the same model with added noise at different levels.
% 
%
% Output:  a row of t-scores as a function of the time shift (from -1 to 1 second)
% 
% noise_vector:  Gaussian noise added to the model to form the data.
%                determined by 'noise_level' .  The noise is multiplied from
%                to give variances from 0.3 to 3 , to build a family of 
%                curves with different noise.
%              
% tau:  width parameter of the HRF used to generate the model of the response
% iti:  mean interscan interval
%
% The function is then shifted in time for ten  different
% delays
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%



%tau = 1.6;
%delta_tau = 0.2;
t_score=[];

%DURATION = 2400; % seconds
DURATION = 600;
SECONDS = 10; % number of points in one second
TR = 1.0;

% create a set of stimuli occurring at random points in the time series
stim=zeros(DURATION*SECONDS,1);
times = (1:iti:DURATION);
times = (times + iti*randn(size(times))) * SECONDS;

% events are not allowed to occur at negative times.
if (~isempty(find(times<1)))
   times = times + abs(min(times)) + 1;
end

times = round(times);
min(times);
stim(times) =1;

%subplot 311, plot(stim);
%load 'fake_iti=4.mat'
%stim = stim(:,1);

% create the HRF using SPM's function
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);
%subplot 312, plot(hrf);

model=conv(stim,hrf);
model = model(1:DURATION*SECONDS);
%subplot 313, plot(model);


% Create the data by convolving the stimulus delta train with the 
% physiological HRF, and then
% add gaussian noise to the response data, and resample every TR.  
% This is the 'sampled data'
   
data = model(1: TR*SECONDS: DURATION*SECONDS);

% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(1*SECONDS,1);...
      model; ...
      zeros(1*SECONDS,1)];


% shift the model by a series of delays
% and compute the correlation coefficient between the data and 
% the model.  


y = data + noise ;
   
for delay= -1: 0.1 : 1
      % shift and resample the model every TR 
      x = model(1*SECONDS + delay *SECONDS+1:TR*SECONDS: size(model)-1*SECONDS +delay*SECONDS);
      
      %hold off
      %plot(x)
      %hold on
      %plot(y, 'g')
      %axis([0 100 0 4])
      %pause
   
      t = my_glm(x,y,[1]);

      t_score= [t_score  t];
end
   

return

