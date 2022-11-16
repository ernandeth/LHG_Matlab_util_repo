function result = shift_visual( datafile, onsets_file, c, DURATION,TR)
%
% function result = shift_visual( datafile, onsets_file, contrast, DURATION, TR)
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



% read in the stimulus onset times and create a spike train
times = load (onsets_file);
times = times/1000;
times = times * TR * SECONDS;

stim=zeros(DURATION*SECONDS,1);
stim(times) =1;



% create the HRF using SPM's function
p = [6 16 1 1 6 0 32 ];

hrf = spm_hrf(1/SECONDS, p);
hrf = hrf/max(hrf);

model=conv(stim,hrf);
model = model(1:DURATION*SECONDS);

% we will use the first regressor to explore what happens when 
% there is an unexpected time shift (here there is only one regressor)
shifting_reg = model(:,1);

% zero pad the regressor, so that it can be shifted in time
shifting_reg = ...
   [zeros(PAD*SECONDS,1); ...
      shifting_reg;...
      zeros(PAD*SECONDS,1)];

% resample the response ever TR
model = model(1: TR*SECONDS: DURATION*SECONDS,  : );


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
	   shifted = shifting_reg( 1 + (PAD + delay)*SECONDS :  TR*SECONDS : (DURATION + PAD +delay)*SECONDS);
	   model(:,1) = shifted;
      
       t = my_glm(model, data, c');
      
       t_score = [t_score t] ;
      
  end      
   
  figure	
  shift = [-1000 * PAD :100:1000 * PAD];
  plot(shift, t_score)
  result = [result;  t_score]; 
   
  % save sim_noise.mat resp_delays model correlations
return
	

