%function result = GA_analysis( conditions, noise, tau,iti, NREGRESSORS, c)
%function result = fake_garavan( noise, tau,iti, NREGRESSORS, c)

noise_var = 0.3;
c = [1 0 0 0]';				% contrast
niterations = 1000;

% Create the Fake Design Matrix (model):

load sim3
conditions = [1 2 3 4];
stimList = myDesign;
ISI = 1.2;

HRF = spm_hrf(ISI);
HRF = HRF/ max(HRF);

clear predMatrix;

for i = 1:size(conditions,2)
   regressor = (stimList == conditions(i));
   regressor = conv(regressor, HRF);
   predMatrix(:,i) = regressor;
end
predMatrix = predMatrix(1:size(stimList),:);      	% eliminate extra values
model = predMatrix;

for i = 1:niterations
	noise = make1overf(300, 1/1.2) * sqrt(noise_var);

	%imagesc(model)
	%colormap(gray)
	result = [];

	% Create the response data by weighting all the regressors by a beta parameter
	% adding all the regressors together
	% and adding noise to the result
	%beta = abs(2.5 * randn(size(model,2), 1));
	beta =[1 1 1 1]';
	data = model * beta * c';
	data = sum(data,2);

	data = data + noise ;
   
   t(i,1) = my_glm(model, data,c);
end
   
hist(t,niterations/50);	

