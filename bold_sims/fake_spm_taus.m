function t_score = fake_spm_taus(err, noise, beta)

%global RSS v var_est t beta_est cov_beta

DURATION = 2400; % seconds
SECONDS = 100; % number of SECONDS in one second
NREGRESSORS = 5;
c = [ 1 0 0 0 0 ];
%iti = 4;
%noise_level = 3;
%tau = 1.5;

% Create the Fake Design Matrix (model_resp):
model_resp = [];
model_error = [];

% create the HRFs

%hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
%hrf_error = make_hrf(2.5*SECONDS,(tau+err)*SECONDS, 20*SECONDS);

% Using the hrfs from SPM ...
p = [6 16 1 1 6 0 32 ];
hrf = spm_hrf( 1/SECONDS, p);
hrf = hrf/max(hrf);

p(1) = p(1) + err;
hrf_error = spm_hrf(1/SECONDS, p);
hrf_error = hrf_error/max(hrf_error);


%stim=[];
%for i=1:NREGRESSORS
  
% create a set of 50 stimuli occurring at random SECONDS in the time series
% the stimuli must be between 8 and 16 seconds apart.
% iti = iti + 1;
   
% s=zeros(DURATION*SECONDS,1);
% times = (1:iti:DURATION-iti);
% times = SECONDS*(times + 4*rand(size(times)));
% s(times) =1;
% stim = [stim s];
%end


% Use a standard set of stimulus onsets instead ...
load 'fake_iti=4.mat'



% convolve the input and the HRF    
% The model response is our prediction (design matrix)
model_resp=conv2(stim,hrf);
model_resp = model_resp(1:DURATION*SECONDS, :);
%whos model_resp stim

model_error = conv2(stim, hrf_error);
model_error = model_error(1:DURATION*SECONDS, :);

% we will use the first regressor to explore what happens when 
% there is an unexpected time shift
shifting_reg = model_error(:,1);
% zero pad the regressor, so that it can be shifted in time
pad = zeros(1*SECONDS, 1);
shifting_reg = ...
   [pad;...
   shifting_reg;...
   pad];


% resample the response every TR
model_resp = model_resp(1: 2*SECONDS : DURATION*SECONDS,  : );
model_error = model_error(1: 2*SECONDS : DURATION*SECONDS,  : );

%whos model_resp model_error noise
%imagesc(model_resp)
%colormap(gray)
%pause

% Create the response data by weighting all the regressors by a beta parameter
% adding all the regressors together
% and adding noise to the result
data = model_resp * beta + noise;
ideal = model_error*beta ;

%plot (data)
%axis([1150 1200 0 30])
%hold on
%plot (ideal,'g')
%pause
%hold off

% shift the regressor by a series of delays: -1000 to 1000 msec.
% and compute the glm parameters for the data and the model
t_score = [];
%figure

for delay= -1: 0.1 : 1

   % shift and resample the the regressor every TR 
   % and put it back in the matrix
   shifted_reg = shifting_reg(1*SECONDS + delay *SECONDS +1: 2*SECONDS : size(shifting_reg)-1*SECONDS + delay*SECONDS);
   model_error(:,1) = shifted_reg;
   
   %whos
   %plot (shifted_reg)
   %axis([0 30 0 3])
   %hold on
   %drawnow
   
   t = my_glm(model_error, data, c');
   t_score = [t_score t] ;  
   
end

%figure
%shift = [-1000:100:1000];
%plot(shift, t_score)
%axis([-500 500 max(t_score)-4 max(t_score) + 0.3 ]);
%drawnow

% save sim_noise.mat resp_delays model_resp correlations
return


