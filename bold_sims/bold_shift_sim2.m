function t_score = bold_shift_sim2(noise, tau, iti, delta_iti)
% function t_score = bold_shift_sim2(noise_vector, tau, iti, delta_iti)
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM ITIs (with a Gaussian distribution), 
% It then computes the t-score of a General Linear Model REgression between the 
% model and the same model with added noise at different levels.
% 
% The program introduces a discrepancy between the model's onset times
% and the data's onset times.  This discrepancy has a gaussian distribution
% with variance determined by delta_iti
%
% Output:  a row of t-scores as a function of the time shift (from -1 to 1 second)
% 			
%
% noise_vector:  Gaussian noise added to the model to form the data.
%                determined by 'noise_level' .  The noise is multiplied from
%                to give variances from 0.3 to 3 , to build a family of 
%                curves with different noise.
%              
% tau:  width parameter of the HRF used to generate the model of the response
% iti:  mean interscan interval
% delta_iti:  variance of discrepancy between the model's ITI and the data's ITI 
%            due to variations in response times...etc.

% The function is then shifted in time for five different
% delays
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%

%tau = 1.6;
%delta_tau = 0.2;

DURATION = 600; % seconds
SECONDS = 10; % number of points in one second
TR = 1;

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

% Create a new set of stimulus times with an additional jitter due to the error
% The data has to be created by convolving the stimulus delta train with the 
% physiological HRF, and then
% add gaussian noise to the response data, and resample every TR.  
% This is the 'sampled data'

if (delta_iti ~=0)
   
   done=0;
   while done==0
      jitters = delta_iti*SECONDS*randn(size(times));
      times = times + jitters;
      
      % again, remember that times are not allowed to be negative ...
      % check to see that is the case ...
      if (isempty(find(times <1)))
         done=1;
      else
         times = times-jitters;
         done=0;
      end
      
   end
   
   times = round(times);
   stim=zeros(DURATION*SECONDS,1);
   stim(times) =1;
  
   data = conv(stim,hrf);   
else
   data = model;
end

data = data(1: TR*SECONDS: DURATION*SECONDS);
y = data + noise ;


% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(1*SECONDS,1);...
      model; ...
      zeros(1*SECONDS,1)];


% shift the model by a series of delays
% and compute the correlation coefficient between the data and 
% the model.  

t_score= [];
for delay= -1: 0.1 : 1
      % shift and resample the model every TR 
      x = model(1*SECONDS + delay *SECONDS+1: TR*SECONDS: size(model)-1*SECONDS +delay*SECONDS);
%hold off
%plot(x)
%hold on
%plot(data, 'g')
%axis([0 100 0 4])
%pause
         
      t = my_glm(x,y,[1]);
      
      t_score= [t_score  t];
      
end




return

