function result = shift_garavan( datafile, regressors, c, DURATION, shifter, drift_file)
%
% function result = shift_garavan( datafile, regressors, contrast, DURATION, shifter, drift_file)
% 
% 

SECONDS = 100; % number of points in one second
PAD = 5;

% create the HRFs
% hrf = make_hrf(2.5*SECONDS,tau*SECONDS, 20*SECONDS);
p = [6 16 1 1 6 0 32 ];
%p(1) = tau;
hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

% Create the Fake Design Matrix (model):
stim = [];
load (regressors)
NREGRESSORS = size(cleanstarts, 2);
stim=zeros(DURATION*SECONDS, NREGRESSORS);
  
for i=1:NREGRESSORS
   stim(cleanstarts(:,i) *2* SECONDS, i) =1;
end

model=conv2(stim,hrf);
model = model(1:DURATION*SECONDS, :);


% we will use the first regressor to explore what happens when 
% there is an unexpected time shift
shifting_reg = model(:,shifter);
% zero pad the regressor, so that it can be shifted in time
shifting_reg = ...
   [zeros(PAD*SECONDS,1); ...
      shifting_reg;...
      zeros(PAD*SECONDS,1)];

% resample the response ever TR
model = model(1: 2*SECONDS: DURATION*SECONDS,  : );

%Append a global mean regressor to the design
if nargin == 6
	drift = load(drift_file);
	drift = drift(:,2) / mean(drift(:,2));
	model = [model  drift];
	NREGRESSORS = NREGRESSORS +1;
	c = [c 0];
end 


imagesc(model)
colormap(gray)
result = [];


% Read the data from file
data = load(datafile);
data = data(:,2);

   
   % shift the regressor by a series of delays: -1000 to 1000 msec.
	% and compute the glm parameters for the data and the model
	t_score = [];
   
   
   for delay= -PAD: 0.1 : PAD
	   
	   % shift and resample the the regressor every TR 
	   % and put it back in the matrix
	   shifted = shifting_reg( 1 + (PAD + delay)*SECONDS :  2*SECONDS : (DURATION + PAD +delay)*SECONDS);
	   model(:,shifter) = shifted;
      t = my_glm(model, data, c');
      t_score = [t_score t] ;
      
      if delay==0
         save my_model.mat model
      end
      
      %whos
      %(1+delay)*SECONDS
      
     % plot(shifted)
     % hold on
     % plot(data,'r')
     % hold off
     % axis([0 50 0 2.5])
  end      
   

  figure
	shift = [-PAD * 1000:100:PAD*1000];
	plot(shift, t_score)
   
  
	result = [result;  t_score]; 
   
   % save sim_noise.mat resp_delays model correlations

   return
	

