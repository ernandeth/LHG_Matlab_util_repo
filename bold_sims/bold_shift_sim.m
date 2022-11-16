function t_score = bold_shift_sim(noise, tau, iti, delta_tau)
% function t_score = bold_shift_sim(noise, tau, iti, delta_tau)
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
% noise_vector:  Gaussian noise added to the model to form the data.
%                determined by 'noise_level' .  The noise is multiplied from
%                to give variances from 0.3 to 3 , to build a family of 
%                curves with different noise.
%              
% tau:  width parameter of the HRF used to generate the model of the response
% iti:  mean interscan interval
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
%subplot 313, plot(model);


% Create the data by convolving the stimulus delta train with the 
% physiological HRF, and then
% add gaussian noise to the response data, and resample every TR.  
% This is the 'sampled data'
if delta_tau ==0
   data = model(1: TR*SECONDS: DURATION*SECONDS);
else
   p(1) = tau + delta_tau;
   p(2) = p(2) + delta_tau;
	hrf = spm_hrf(1/SECONDS, p);
	hrf = hrf/max(hrf);
   data = conv(stim,hrf);   
   data = data(1: TR*SECONDS: DURATION*SECONDS);
end

% zero pad the model to allow the correlation to be shifted in time
% adding 1 second worth of samples
model = ...
      [zeros(1*SECONDS,1);...
      model; ...
      zeros(1*SECONDS,1)];


% shift the model by a series of delays
% and compute the correlation coefficient between the data and 
% the model.  

% DO this 10 times with 10 different noise levels

% rhos_matrix = [];
%t_score_matrix=[];

%for i = 0.3 : 0.3 :3
   % recall that variance will increase as the square root...
   
    %y = data + sqrt(i)*noise ;

    y = data + noise ;
   
   t_score= [];
   rhos = [];
   for delay= -1: 0.1 : 1
      % shift and resample the model every TR 
      x = model(1*SECONDS + delay *SECONDS+1: TR*SECONDS: size(model)-1*SECONDS +delay*SECONDS);
      
      %hold off
      %plot(x)
      %hold on
      %plot(y, 'g')
      %axis([0 100 0 4])
      %pause
   
      t = my_glm(x,y,[1]);
      %t = my_glm(x,y,[1]);
      %rho = corrcoef(x,y);

      %rhos = [rhos rho(1,2)];
      t_score= [t_score  t];
   end
   
    %rhos_matrix = [rhos_matrix;  rhos];
    %t_score_matrix = [t_score_matrix ; t_score];
      
    %end


% convert correlation coefficent to t score:
%df = size(data,1) - 2;
%t_score = rhos_matrix .* sqrt(df) ./ sqrt(1 - rhos_matrix.^2);

%t_score = t_score_matrix;

return

