function result = fake_garavan2( tau,data_file, c)
%function result = fake_garavan2(tau,data_file, c)

%DURATION = 2400; % seconds
DURATION = 968*2; % seconds
SECONDS = 10; % number of points in one second
TR = 2;

% create the HRFs
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

% Create the Fake Design Matrix (model):
stim = [];
load regressors7
NREGRESSORS = size(cleanstarts, 2);
stim=zeros(DURATION*SECONDS,NREGRESSORS);
  
for i=1:NREGRESSORS
   stim(cleanstarts(:,i) *TR* SECONDS, i) =1;
end

model=conv2(stim,hrf);
model = model(1:DURATION*SECONDS, :);

% we will use the first regressor to explore what happens when 
% there is an unexpected time shift
shifting_reg = model(:,4);
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

% Read the data from file
data = load(data_file);
data = data(:,2);


   % shift the regressor by a series of delays: -1000 to 1000 msec.
	% and compute the glm parameters for the data and the model
	t_score = [];
   
   
   for delay= -1: 0.1 : 1
	   
	   % shift and resample the the regressor every TR 
	   % and put it back in the matrix
	   shifted = shifting_reg( 1 + (1 + delay)*SECONDS :  TR*SECONDS : (DURATION + 1 +delay)*SECONDS);
	   model(:,4) = shifted;
      t = my_glm(model, data,c');
      t_score = [t_score t] ;
      
      %whos
      %(1+delay)*SECONDS
      
      if ((abs(delay) < 0.001) & 1)
          plot(x*100)
          hold on
          plot(data - mean(data), 'g')
          %axis([0 200 0 4])
          drawnow
          hold off
          %pause
       end      
  end      
   
%   figure
%	shift = [-1000:100:1000];
%	plot(shift, t_score)
   
   %axis([-500 500 max(t_score)-4 max(t_score) + 0.3 ]);
   %drawnow
   
	result = [result;  t_score]; 
   
   % save sim_noise.mat resp_delays model correlations
return
	

