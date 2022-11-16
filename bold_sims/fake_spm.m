function result = fake_spm( noise, tau,iti, NREGRESSORS, c)
%function result = fake_spm( noise, tau,iti, NREGRESSORS, c)

DURATION = 600; % seconds
SECONDS = 10; % number of points in one second
TR = 1;

% create the HRFs
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

% Create the Fake Design Matrix (model):
stim = [];

for i=1:NREGRESSORS
   
  s=zeros(DURATION*SECONDS,1);
  times = (1:iti:DURATION);
  times = (times + (iti/2)*rand(size(times))) * SECONDS;
  s(times) =1;
  stim = [stim s];
   
end

model=conv2(stim,hrf);
model = model(1:DURATION*SECONDS, :);

% we will use the first regressor to explore what happens when 
% there is an unexpected time shift
shifting_reg = model(:,1);
% zero pad the regressor, so that it can be shifted in time
shifting_reg = ...
   [zeros(1*SECONDS,1); ...
      shifting_reg;...
      zeros(1*SECONDS,1)];

% resample the response ever TR
model = model(1: TR*SECONDS: DURATION*SECONDS,  : );


%imagesc(model)
%colormap(gray)
result = [];

% Create the response data by weighting all the regressors by a beta parameter
% adding all the regressors together
% and adding noise to the result
%beta = abs(2.5 * randn(size(model,2), 1));
beta =1;
data = model * beta;

for i= 0.5 : 0.5 : 5
   experiment = data + sqrt(i)*noise ;
   
   % shift the regressor by a series of delays: -1000 to 1000 msec.
	% and compute the glm parameters for the data and the model
	t_score = [];
   
   
   for delay= -1: 0.1 : 1
	   
	   % shift and resample the the regressor every TR 
	   % and put it back in the matrix
	   shifted = shifting_reg( 1 + (1 + delay)*SECONDS :  TR*SECONDS : (DURATION + 1 +delay)*SECONDS);
	   model(:,1) = shifted;
      t = my_glm(model, experiment,c');
      t_score = [t_score t] ;
      %whos
      %(1+delay)*SECONDS
      
     % plot(shifted)
     % hold on
     % plot(data,'r')
     % hold off
     % axis([0 50 0 2.5])
     
%      pause

   end
      
	%figure
	%shift = [-1000:100:1000];
	%plot(shift, t_score)
	%axis([-500 500 max(t_score)-4 max(t_score) + 0.3 ]);
   %drawnow
   
	result = [result;  t_score]; 
   %whos result t_score
   
end

   
	% save sim_noise.mat resp_delays model correlations
return
	

