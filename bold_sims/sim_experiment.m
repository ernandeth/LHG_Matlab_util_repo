function result = sim_experiment( regressors, c, beta, DURATION, TR, noise)
%
%function result = sim_experiment( regressors, contrast, beta, DURATION, TR, noise)
% 
% 

SECONDS = 100; % number of points in one second
beta = beta';

% create the HRFs
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
%p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

% Load the file with the stimulus onset times and 
% Create the Fake Design Matrix (model) by convolving the 
% train of spikes with the HRF
stim = [];
load (regressors)
NREGRESSORS = size(stim_times, 2);
stim=zeros(DURATION*SECONDS, NREGRESSORS);
  
for i=1:NREGRESSORS
   stim(stim_times(:,i) *TR* SECONDS, i) =1;
end

model=conv2(stim,hrf);
model = model(1:DURATION*SECONDS, :);



% resample the response ever TR
model = model(1: TR*SECONDS: DURATION*SECONDS,  : );
figure
imagesc(model)
colormap(gray)
title(' This is the design matrix')
xlabel('regressors or effects')
ylabel('scan number')


result = [];
% Generate some fake data according to the specified matrix and 
% the weights of each regressors, ie- the betas
noise = noise(1:size(model,1));
data = model*beta + noise;

   
% Compute the general linear model stuff
t = my_glm(model, data, c');

% check your results
% Note: recall that my_glm adds an extra regressor for
% DC offset or baseline, so we need to take that off
load glm
beta_est = beta_est(2:end);
data_est = model * beta_est;

figure
plot(data)
hold
plot(data_est, 'r')

t
beta
beta_est

return
	

