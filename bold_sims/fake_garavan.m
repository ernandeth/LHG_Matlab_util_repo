function result = fake_garavan( noise,  NREGRESSORS, c, beta, shifter)
%function result = fake_garavan( noise, NREGRESSORS, c, beta, shifter)

%DURATION = 2400; % seconds
DURATION = 968*2; % seconds
SECONDS = 100; % number of points in one second
tau = 6;
iti = 16;
PAD = 2;


noise = noise(1:DURATION/2);

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


%imagesc(model)
%colormap(gray)
result = [];

% Create the response data by weighting all the regressors by a beta parameter
% adding all the regressors together
% and adding noise to the result
data = model * beta' ;
data = sum(data,2);

   data = data + noise ;
   
   % shift the regressor by a series of delays: -1000 to 1000 msec.
	% and compute the glm parameters for the data and the model
	t_score = [];
   
   
   for delay= -PAD: 0.1 : PAD
	   
	   % shift and resample the the regressor every TR 
	   % and put it back in the matrix
	   shifted = shifting_reg( 1 + (PAD + delay)*SECONDS :  2*SECONDS : (DURATION + PAD +delay)*SECONDS);
	   model(:,shifter) = shifted;
       
      t = my_glm(model, data,c');
      t_score = [t_score t] ;
      %whos
      %(1+delay)*SECONDS
      
     % plot(shifted)
     % hold on
     % plot(data,'r')
     % hold off
     % axis([0 50 0 2.5])
  end      
   
   figure
	shift = [-1000*PAD:100:1000*PAD];
	plot(shift, t_score)
   
   %axis([-500 500 max(t_score)-4 max(t_score) + 0.3 ]);
   %drawnow
   
	result = [result;  t_score]; 
   
   % save sim_noise.mat resp_delays model correlations
return
	

