function t_score = bold_shift_data3(data_file, onsets_file,tau, n_iti)
% function t_score = bold_shift_data3(data_file,onsets_file,tau, n_iti)
% Luis Hernandez
% University of Michigan
%
% This program generates a simulated BOLD response
% to a set of EVENTS that occur AT RANDOM ITIs (with a Gaussian distribution), 
% It then computes the t-score of a General Linear Model REgression between the 
% model and the same model with added noise.
%
% The program introduces discrepancy between the number of events in the 
% model and the data, given by the variable n_iti.
% 
% Output:  a row of t-scores as a function of the time shift (from -1 to 1 second)
% 			
%
% noise_vector:  Gaussian noise added to the model to form the data.
%                determined by 'noise_level' .  The noise is multiplied from
%                to give variances from 0.3 to 3 , to build a family of curves with different noise.             
% tau:  width parameter of the HRF used to generate the model of the response
% iti:  mean interscan interval
% n_iti:  number of events that are missing from the data (no response) or that 
% 			are included in the data by error, as a response to some extraneous process
%
% The function is then shifted in time for five different
% delays
%
% The total imaging time is DURATION sec.
%
% The number of points in one second is specified by the 
% variable 'SECONDS'.( not related to TR )
%

%tau = 6;
%iti=16;
%n_iti = -10;

DURATION = 600; % seconds
SECONDS = 10; % number of points in one second

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

model=conv(stim,hrf);
m = model;

% Read the data from file
data = load(data_file);
data = data(:,2);

% remove n_iti events from the set of stimuli and re-create the data
% The data has to be created by convolving the stimulus delta train with the 
% physiological HRF, and then
% add gaussian noise to the response data, and resample every TR.  
% This is the 'sampled data'


if (n_iti > 0)
   
   % add responses to the model
   times2= size(stim,1) * rand(abs(n_iti),1);
   times2 = ceil(times2);
   stim(times2) = 1;
   
   model = conv(stim,hrf); 
   fprintf('\rAdding Responses ....');

   
elseif (n_iti < 0)
   % remove responses from the model
   missing= 1 + (size(times,1)-1) * rand(abs(n_iti),1);
   missing=floor(missing);

   stim(times(missing)) = 0;
   fprintf('\rRemoving Responses ....');

   model = conv(stim,hrf);   

else
   
   fprintf('\rRight number of Responses ....');
      
end

% Remove the tail end of the convolution
model = model(1:DURATION*SECONDS);
%model = model(1: 2*SECONDS: DURATION*SECONDS);


y = data ;

% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(1*SECONDS,1);...
      model; ...
      zeros(1*SECONDS,1)];

m = m(1:1*SECONDS:DURATION*SECONDS);

% shift the model by a series of delays
% and compute the correlation coefficient between the data and 
% the model.  
t_score= [];

for delay= -1: 0.1 : 1
      % shift and resample the model every TR 
      
      x = model(...
         1*SECONDS + delay*SECONDS + 1: ...
         1*SECONDS: ...
         size(model)- 1*SECONDS + delay*SECONDS);
      
      % some plotting stuff for debugging purposes ... 
           
      if ((abs(delay) < 0.001) & 0)
         %hold off
         % plot(data - mean(data), 'g')
         %hold on
         %whos x m
         plot(x - m)
         %axis([0 200 0 4])
         drawnow
         pause
		end      
      
      %whos x y
      t = my_glm(x,y,[1]);
      
%  rho = corrcoef(x,y);
%  rho = rho(2,1);    
%  df = 1198;
%  t = rho * sqrt(df) ./ sqrt(1 - rho.^2);

      t_score= [t_score  t];
      
end




return

