% trainNetwork_VSI_20210929
fprintf('\n Loading training data ...');
load chirpfun_vss/trainData.mat
load chirpfun_vss/trainTgts.mat

% clip the first 10 frames - they are not very good:
% trainData = trainData(:, 11:end);

% Normalize training data so that the first value is always 1
fprintf('\n Scaling training data ...');
for n=1:size(trainData,1)
%    trainData(n,:) = trainData(n,:) / trainData(n,1);
%    trainData(n,:) = trainData(n,:) / trainData(n,end);

    % 11/22/21 : try normalizing the training data, instead of scaling to
    % the first sample
    trainData(n,:) = trainData(n,:) / norm(trainData(n,:));
end

% Add noise - (warning - training data already had noise)
sigma = 0.0025;
sigma = 0.001;

noiseMat = randn(size(trainData)) * sigma;
trainData = trainData + noiseMat;

Ntrials = size(trainData,1);
Nframes = size(trainData,2);

% order of the targets (params to estimate)
parm_number = 1;

parm_names = {
    'cbf';
    'cbv';
    'kfor'; % not used
    'bat2';  % not used
    'bat';
    'r1';
    'flip';
    'r2tiss';
    };  % these correspond to the 'targets' matrix columns

% scale targets so that they average to 1 - keep track of the scaling !!
fprintf('\n Scaling training parameter (targets) %d ...', parm_number);
parm_train = trainTgts(:,parm_number);
scale = 1/mean(parm_train);  % scale so that the average of the population is 1
parm_train = parm_train * scale;

% get the training data in the right shape for the network
trainData = reshape(trainData', Nframes, 1,1,Ntrials);
%%
%{
% first prep the GPU
% compiling the libraries - this may be needed only the first time:
fprintf('\n Setting up GPU ...');
parallel.gpu.enableCUDAForwardCompatibility(true);
mygpu = gpuDevice(1);
%}
%%

% create network 1
fprintf('\n Creating Network ...');
%{
trainData = squeeze(trainData);
Mynet = [
        sequenceInputLayer(Nframes);
        lstmLayer(30, 'OutputMode', 'Last');
        fullyConnectedLayer(1);
        regressionLayer;
        ];
%}
%
Mynet = [
    imageInputLayer([Nframes 1]);
    % convolution2dLayer([3 1], 2,'stride',1);  % no good for us
    
    fullyConnectedLayer(floor(Nframes));
    reluLayer();
    
    fullyConnectedLayer(floor(Nframes));
    reluLayer();
    
    fullyConnectedLayer(floor(Nframes));
    reluLayer();
        
    fullyConnectedLayer(floor(Nframes/2));
    reluLayer();
    
    fullyConnectedLayer(floor(Nframes/4));
    reluLayer();
    
    fullyConnectedLayer(floor(Nframes/4));
    reluLayer();
    
   fullyConnectedLayer(floor(Nframes/8));
   reluLayer();
    
   fullyConnectedLayer(floor(Nframes/8));
   reluLayer();
    
    fullyConnectedLayer(1);
    regressionLayer;
    ]
%}
%
% train Network
%
opts = trainingOptions('adam', ... % also try sgdm
    'MaxEpochs', 20, ...  % 20
    'Plots', 'training-progress', ...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-5, ... %1e-4
    'WorkerLoad', ones([24 1])) ;
% things to try with ADAM
%'GradientDeacayFactor', 0.85,
%'SquaredGradientDeacayFactor', 0.85,

%    'OutputNetwork', 'best-validation-loss', ... % this is for sgdm

fprintf('\n Training Network ...');
[Mynet diagnostic]  = trainNetwork(...
    (trainData) , parm_train, Mynet, opts );

str=sprintf('save Mynet_%s_chirpfun_vss.mat Mynet scale', parm_names{parm_number})
eval(str)
%
% testing
%
fprintf('\n Loading Test Data ...');

load chirpfun_vss/testData.mat
load chirpfun_vss/testTgts.mat

fprintf('\n loading network Mynet_%s_chirpfun_vss.mat ', parm_names{parm_number})
str=sprintf('load Mynet_%s_chirpfun_vss.mat ', parm_names{parm_number})
eval(str)

% clip the first 10 frames - they are not very good:
% testData = testData(:, 11:end);

Ntrials = size(testData,1);
Nframes = size(testData,2);

fprintf('\n Scaling Test Data ...');
% Normalize testing data so that the first value is always 1
for n=1:size(testData,1)
    %testData(n,:) = testData(n,:) / testData(n,1);
    %testData(n,:) = testData(n,:) / testData(n,end);
    
    % 11/22/21 : try normalizing the training data, instead of scaling to
    % the first sample
    testData(n,:) = testData(n,:) / norm(testData(n,:));
end

testData = reshape(testData', Nframes, 1,1,Ntrials);
est = zeros(size(testTgts,1),1);


fprintf('\n Testing the Network ...');

parfor n=1:Ntrials

    est(n) = predict(Mynet, testData(:,:,:, n));
        
    if mod(n,1e4)==0
        fprintf('\nTesting time course number ... %d (of %d)', n, Ntrials);
    end
end

% use the target scaling factor to scale the outputs into the apporpriate range
est = est/scale;
true = testTgts(:,parm_number); 

% now compare estimates versus the truth 
figure
plot(true, est, '.'); 
ylabel('Prediction')
xlabel('Truth');
title(parm_names{parm_number})

axis([-0.01 2 -0.01 2]/scale)
rho = corrcoef(est, true);

% fit quadratic model.  linearity? bias?  high order?
[p s] = polyfit(true, est,2);
x= linspace(0, max(true), 100);
hold on
plot(x, p(1)*x.^2 + x*p(2) + p(3), 'g')
hold off

legend (sprintf('correlation : %02f', rho(1,2)), ...
    sprintf('%02f * x^2 + %02f * x  + %02f',p), '');
hold on; plot([-100 100], [-100 100]) ; hold off




