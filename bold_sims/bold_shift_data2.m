function t_score = bold_shift_data2(data_file, onsets_file,tau, delta_iti)
% function t_score = bold_shift_data2(data_file,onsets_file,tau, delta_iti)
%
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM ITIs (with a Gaussian distribution), 
% It then computes the t-score of a General Linear Model REgression between the 
% model and the same model with added noise at different levels.
% 
% The program introduces a discrepancy between the time constants of the 
% model's underlying HRF, and those of the data's HRF. the discrepancy is 
% specified by delta_tau
%
% Output:  a row of t-scores as a function of the time shift (from -1 to 1 second)
% 			
% data_file: contains a timeseries from a region 
%           (2 columns: number intensity)           
% tau:  width parameter of the HRF used to generate the model of the 
%       response
% onsets_file:  contains the onset times of the experiment
% delta_iti:  variance of discrepancy between the model's ITI and the data's ITI 
%            due to variations in response times...etc.
%
% The function is then shifted in time for five different
% delays
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%


DURATION = 600; % seconds
SECONDS = 10; % number of points in one second

% read in the stimulus onset times and create a spike train
times = load (onsets_file);
times = round(times / 1000);
times = times(find(times < DURATION)) * SECONDS;

stim=zeros(DURATION*SECONDS,1);
stim(times) =1;
%subplot 311, plot(stim);
%load 'fake_iti=4.mat'
%stim = stim(:,1);

% create the HRF using SPM's function
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau ;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

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
end



model=conv(stim,hrf);
model = model(1:DURATION*SECONDS);

%subplot 313, plot(model);

% Read the data from file
data = load(data_file);
data = data(:,2);

% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(1*SECONDS,1);...
      model; ...
      zeros(1*SECONDS,1)];


% shift the model by a series of delays
% and compute the correlation between the data and 
% the model.  

y = data;
t_score= [];
   
for delay= -1: 0.1 : 1
   
   % shift and resample the model every TR 
   x = model(1*SECONDS + delay*SECONDS + 1 : ...
      1*SECONDS : ...
      size(model)-1*SECONDS + delay*SECONDS);
   
   if ((abs(delay) < 0.001) & 0)
      
      hold off
      plot(x)
      hold on
      plot(data - mean(data), 'g')
      %axis([0 200 0 4])
      drawnow
      pause
      
   end      
   
   t = my_glm(x,y,[1]);
   t_score= [t_score  t];
   
end

return

