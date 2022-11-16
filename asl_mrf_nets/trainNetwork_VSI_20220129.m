%%
%{
% first prep the GPU
% compiling the libraries - this may be needed only the first time:
fprintf('\n Setting up GPU ...');
parallel.gpu.enableCUDAForwardCompatibility(true);
mygpu = gpuDevice(1);
%}
%%

% trainNetwork_VSI_20210929
fprintf('\n Loading training data ...');
load AL_schedule_model_220130/trainData_3.mat
load AL_schedule_model_220130/trainTgts_3.mat

% clip the first 10 frames - they are not very good:
%trainData = trainData(:, end-190:end); 2/17/22
trainData = trainData(:, 2:200);
%trainData = trainData(:, 3:200); % 2/21/22


% Normalize training data so that the first value is always 1
fprintf('\n Scaling training data ...');
for n=1:size(trainData,1)
    trainData(n,:) = trainData(n,:) / mean(abs(trainData(n,:)));
end

Ntrials = size(trainData,1);
Nframes = size(trainData,2);

% Add noise - (warning - training data already had noise)
sigma = 0.05;
%sigma = 1e-3; % 2/19/22
%sigma = 0;


noiseMat = randn(size(trainData)) ;

%2/21/22 : try adding AR1 noise
%ARmat = makeAR1mat(0.5, Nframes);
%ARmat = makeAR1mat(0.9, Nframes);
ARmat = makeAR1mat(0.7, Nframes);

noiseMat = ARmat * noiseMat';
noiseMat = noiseMat';

noiseMat = noiseMat/norm(noiseMat);

trainData = trainData + noiseMat*sigma;

% 2/28/22 - Data augmentation with pairwise difference to see if it helps
%{
fprintf('\n Extending training data ...');
myDiffs = trainData(:, 2:2:end) - trainData(:, 1:2:end-1);
trainData = [trainData myDiffs];
%}

Nframes = size(trainData,2);
% get the training data in the right shape for the network
trainData = reshape(trainData', Nframes, 1,1,Ntrials);

% Now get the TESTING DATA and do the same scaling etc
fprintf('\n Loading Test Data ...');

load AL_schedule_model_220130/testData_3.mat
load AL_schedule_model_220130/testTgts_3.mat

% clip the first 10 frames - they are not very good:
%testData = testData(:, end-190:end);
testData = testData(:, 2:200); 
%testData = testData(:, 3:200);  % 2/21/22

% 2/19/22 - do the difference to see if it helps
%testData = diff(testData,1,2);

% 2/21/22 - pairwise difference to see if it helps
%testData = testData(:, 2:2:end) - testData(:, 1:2:end);

Ntests = size(testData,1);
Nframes = size(testData,2);

fprintf('\n Scaling Test Data ...');
% Scale  testing data so that the first value is always 1
for n=1:size(testData,1)
    testData(n,:) = testData(n,:) / mean(abs(testData(n,:)));
end
%}

%2/21/22 : try adding AR1 noise
noiseMat = randn(size(testData)) ;

noiseMat = ARmat * noiseMat';
noiseMat = noiseMat';

noiseMat = noiseMat/norm(noiseMat);

testData = testData + noiseMat*sigma;

% 2/28/22 - Data augmentation with pairwise difference to see if it helps
%{
fprintf('\n Extending testing data ...');
myDiffs = testData(:, 2:2:end) - testData(:, 1:2:end-1);
testData = [testData myDiffs];
%}

Nframes = size(testData,2);
testData = reshape(testData', Nframes, 1,1,Ntests);

whos train* test*

%% order of the targets (params to estimate)
parm_names = {
    'cbf';
    'cbv';
    'bat';
    'r1';
    'r2';
    'flip'
    };  % these correspond to the 'targets' matrix columns
%%

for parm_number = 1:6
    
    % scale targets so that they average to 1 - keep track of the scaling !!
    fprintf('\n Scaling training parameter (targets) %d ...', parm_number);
    parm_train = trainTgts(:,parm_number);
    scale = 1/mean(parm_train);  % scale so that the average of the population is 1
    parm_train = parm_train * scale;
    
    Mynet = [];
    % create network 
    fprintf('\n Creating Network ...');
    Mynet = [
        imageInputLayer([Nframes 1]);
        %convolution2dLayer([3 1], 2,'stride',1);
        
        fullyConnectedLayer(floor(Nframes));
        reluLayer();
        
        fullyConnectedLayer(floor(Nframes));
        reluLayer();
       
        fullyConnectedLayer(floor(Nframes)); % 2/19/22
        reluLayer();

        fullyConnectedLayer(floor(Nframes)); % 2/23/22
        reluLayer();
        
%         fullyConnectedLayer(floor(Nframes)); % 3/17/22 for CBV
%         reluLayer();
%         
%         fullyConnectedLayer(floor(Nframes/4));
%         reluLayer();
%         
%         fullyConnectedLayer(floor(Nframes/4));
%         reluLayer();
%         
%         fullyConnectedLayer(floor(Nframes/8));
%         reluLayer();
%         
%         fullyConnectedLayer(floor(Nframes/8));
%         reluLayer();
        
        fullyConnectedLayer(1);
        regressionLayer;
        ]
    
    %
    % train Network
    %
     
    Nepochs = 10;
    if parm_number==3
        Nepochs=20;
    end
   
    opts = trainingOptions('adam', ... % also try sgdm
        'MaxEpochs', Nepochs, ...
        'Shuffle','every-epoch', ...
        'InitialLearnRate', 1e-5, ... %1e-4
        'ValidationPatience', 100, ...
        'WorkerLoad', ones([24 1])) ;
    %'OutputNetwork', 'best-validation-loss', ... 
    %  'Plots', 'training-progress', ...

    fprintf('\n Training Network ...');
    [Mynet diagnostic]  = trainNetwork(...
        trainData , parm_train, Mynet, opts );
    %
    str=sprintf('save AL_schedule_model_220130/Mynet_%s.mat Mynet scale', parm_names{parm_number})
    eval(str)
    %
    % testing
    %
    
    fprintf('\n loading network Mynet_%s.mat ', parm_names{parm_number})
    str=sprintf('load  AL_schedule_model_220130/Mynet_%s.mat ', parm_names{parm_number})
    eval(str)
    
    est = zeros(size(testTgts,1),1);
    
    
    fprintf('\n Testing the Network ...');
    
    parfor n=1:Ntests
        
        est(n) = predict(Mynet, testData(:,:,:, n), 'ExecutionEnvironment','cpu');
        
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
    
    drawnow
end

