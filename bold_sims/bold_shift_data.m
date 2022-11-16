function t_score = bold_shift_data(data_file, onsets_file,tau, delta_tau, DURATION, TR)
% function t_score = bold_shift_data(data_file,onsets_file,tau, delta_tau, DURATION, TR)
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
% delta_tau:  Error in the Tau parameter for the Model.  The data are 
%           generated 
%           using tau = tau+ delta_tau
%
% The function is then shifted in time for five different
% delays
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%


%DURATION = 600; % seconds
%DURATION = 968*2; % seconds
SECONDS = 10; % number of points in one second

% read in the stimulus onset times and create a spike train
times = load (onsets_file);
%times = round(times / 1000);
%times = times(find(times < DURATION)) * 2* SECONDS;
times = times * TR * SECONDS;

stim=zeros(DURATION*SECONDS,1);
stim(times) =1;
%subplot 311, plot(stim);
%load 'fake_iti=4.mat'
%stim = stim(:,1);

% create the HRF using SPM's function
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau + delta_tau;
p(2) = p(2) + delta_tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

model=conv(stim,hrf);
model = model(1:DURATION*SECONDS);
%subplot 313, plot(model);

% Read the data from file
data = load(data_file);
data = data(:,2);

% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(10*SECONDS,1);...
      model; ...
      zeros(10*SECONDS,1)];


% shift the model by a series of delays
% and compute the correlation between the data and 
% the model.  

y = data;
t_score= [];
   
for delay= -10: 0.1 : 10
   
   % shift and resample the model every TR 
   x = model(10*SECONDS + delay*SECONDS + 1 : ...
      TR*SECONDS : ...
      size(model)-10*SECONDS + delay*SECONDS);
   
   if ((abs(delay) < 0.001) & 1)
      plot(x*100)
      hold on
      plot(data - mean(data), 'g')
      %axis([0 200 0 4])
      drawnow
      hold off
      %pause

   end      

   t = my_glm(x,y,[1]);
   t_score= [t_score  t];
   
end

return

