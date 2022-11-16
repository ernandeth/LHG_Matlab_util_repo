clear all
%
addpath  /home/hernan/matlab/flow/fingerprint
read_parms_file = 0

dt = 1e-3; % <----  note 


% create ASL sequence timing
if read_parms_file
    % read in schedule from file:
    aq_parms.t_tag          = load('t_tags.txt');
    aq_parms.t_delays       = load('t_delays.txt');
    aq_parms.t_adjusts      = load('t_adjusts.txt');;
    aq_parms.labelcontrol   = load('isVelocitySelective.txt');
    aq_parms.doArtSup       = load('doArtSuppression.txt');;
    aq_parms.ArtSup_delay   = 0.1 *ones(size(delays)); % delay between AS pulse and acqusition
    aq_parms.t_aq           = load('t_aqs.txt');
    Nframes = length(t_delay);
else
    % alternative schedule:
    %
    Nframes = 100;
    idx = randperm(Nframes);
    
    aq_parms.labelcontrol           = ones(Nframes,1);
    aq_parms.labelcontrol(2:2:end)  = 0;
    % aq_parms.labelcontrol(idx(1:end/2))  = 0;    % randomize labels
    % aq_parms.labelcontrol(1:8)        = -1;
     
    
    aq_parms.doArtSup               = ones(Nframes,1);
    % aq_parms.doArtSup(5:5:end/2)    = 0;
    % aq_parms.doArtSup(1:4)          = -1;
    
    aq_parms.t_aq                   = 0.5 *ones(Nframes,1);
    aq_parms.ArtSup_delay           = 0.1 *ones(Nframes,1); % delay between AS pulse and acqusition
    aq_parms.t_adjusts              = 1.5 *ones(Nframes,1);
    aq_parms.t_tag                  = zeros(Nframes,1);
    %aq_parms.t_delays               =  0.6 + 0.9*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    aq_parms.t_delays               =  1.2 + 0.3*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    aq_parms.t_delays               =  0.1 + 1.1*abs(chirp(linspace(0,1,Nframes), 0, 1,4))';
    %aq_parms.t_delays(2:2:end)      = aq_parms.t_delays(1:2:end);


end

% default parameters
parms.f         = 60/6000;
parms.cbva      = 0.02;
parms.bat       = 0.5;
parms.r1tis     = 1/1.4;

% fixed defaults
parms.flip      = deg2rad(90);
parms.alpha_ai  = 0.8; % arterial inversion efficiency
parms.alpha_ti  = 0.75; % tissue inversion efficiency
parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression
parms.Mtis0     = 1;

%
% sanity check: show the signal for one case
doSub = 0;
doFigs = 0;
figure
subplot(211)
parms.f = 60/6000;
parms.bat = 0.5;
%entry1 = abs(gen_signals_vs_201020(parms, aq_parms, doFigs,doSub, dt));
entry1 = abs(gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
plot(entry1/entry1(1)); 
hold on

% parms.r1tis = 0.7*parms.r1tis;
%parms.f = 30/6000;
parms.bat = 0.9;

%entry2 = abs(gen_signals_vs_201020(parms, aq_parms, doFigs,doSub, dt));
entry2 = abs(gen_signals_vs_201215(parms, aq_parms, doFigs,doSub, dt));
plot(entry2/entry2(1)); 
legend('f = 60', 'f = 30')
legend('bat = 0.5', 'bat = 0.9')
hold off

subplot(212)
plot((entry1-entry2)./entry1);
legend('subtraction')
drawnow


%%

fprintf('\nGenerating training Parameters ....');
% scale the training parameters so that the weights can be adjusted in the
% RMSE calculation at the end of the network  during backpropagation
n = 1;
Ntrain =  5e6; %1.2e7;      % memory used for training parms: Ntrain x 7 x 4 Bytes.

% memory for training data is : Ntrain x Nframes x 4 Bytes
% This means if Ntrain = 10^7
% training parms size is 280 MB
% training data size is 12,000 MB

% again, for Anish's test:
Nparms = 6;

fvals =    single( rand(Ntrain, 1) * 0.02);     % ml/ml/s

cbvavals = single( rand(Ntrain, 1) * 0.05 );     % fraction
%cbvavals(:) = 0.01;
batvals = round(100*single( (rand(Ntrain, 1) * 0.1 + 0.9) ))/100;     % s
%batvals(:) = 0.5;
r1vals =   single( (rand(Ntrain, 1)* 2 + 0.3));     % 1/s
%r1vals(:) = 1/1.4;

flipvals = single( pi/2*(0.995 + 0.01*(rand(Ntrain, 1))));     % rad : 90 +/- 10% 
flipvals(:) = pi/2;

alpha_tivals  = single(rand(Ntrain, 1)*0.1 + 0.7 );  % tissue inversion efficiency  <--- reduced variability.
alpha_tivals(:) = 0.8;

% allocate space and store training parms into single matrix
Nframes = length(entry1);

Ncombinations = Ntrain;

trainingData = single(zeros(Ncombinations, Nframes));

%trainingParms = [fvals cbvavals batvals  r1vals flipvals alpha_tivals];
trainingParms = [fvals cbvavals batvals  r1vals flipvals alpha_tivals];


%  synth. training data
fprintf('\nGenerating training data ....');
tic
doFigs = 0;
doSub = 0;
parfor n=1:Ncombinations
    
    % fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    tparms = parms;

   
    tparms.f =       trainingParms(n, 1) ;
    tparms.cbva =    trainingParms(n, 2) ;
    tparms.bat =     trainingParms(n, 3) ;
    tparms.r1tis =   trainingParms(n, 4) ;
    
    % exclude these from training for now: use defaults
    % parms.flip =    trainingParms(n, 5) ;
    % parms.alpha_ti= trainingParms(n, 6) ; % tissue inversion efficiency
    %
    entry = abs(gen_signals_vs_201215(tparms, aq_parms, doFigs, doSub, dt));
    %entry = abs(gen_signals_vs_201020(tparms, aq_parms, doFigs, doSub, dt));
    entry = entry/entry(1);
    %plot(entry);drawnow
    
    trainingData(n,:) = entry;
    
end
toc
trainingData = reshape(trainingData', Nframes, 1,1,Ncombinations);
% add noise to the training data
% trainingData = trainingData + randn(size(trainingData))*0.01;
whos

%%
% try differencing the data:
% D_trainingData = trainingData(2:2:end, :,:,:) - trainingData( 1:2:end, :,:,:);
% Nframes = Nframes/2;

%
fprintf('\nGenerating Network ....');

inputLayer = imageInputLayer([Nframes 1]);

%bn =batchNormalizationLayer;

Mynet = [
    imageInputLayer([Nframes 1]);
    convolution2dLayer([3 1], 2,'stride',1);
    fullyConnectedLayer(Nframes);
    reluLayer();
    fullyConnectedLayer(Nframes/2);
    reluLayer();
%     fullyConnectedLayer(round(Nframes/4));
%     reluLayer();  
    fullyConnectedLayer(1);
    regressionLayer;
    ]

% Try including batch normalization:
% Mynet=[inputLayer;   f1; r1; bn;   f2; r1; bn;   f3; r1; bn; f4;  r6]

%delete(gcp)

opts = trainingOptions('adam', ... % do not use adam for perfusion
    'MaxEpochs', 10, ...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-2, ...
    'WorkerLoad', ones([40 1]) )
        
    % 'Shuffle','every-epoch', ...
    %'Plots', 'training-progress', ...
    %'GradientDecayFactor', 0.19000 , ...
    %'SquaredGradientDecayFactor', 0.19990 , ...   'ExecutionEnvironment', 'parallel', ...
    
    %'ValidationPatience', 5, ...
    %'MiniBatchSize', 500, ...
    %'InitialLearnRate', 1e-4, ...
    %'LearnRateSchedule', 'piecewise', ...


save NN_stuff Mynet aq_parms Nframes  trainingData trainingParms

noise_level =  0; %1e-7;

%{
fprintf('\nTraining CBF Network ....');

[Mynet_cbf diagnostic] = trainNetwork(...
    trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,1), ...
    Mynet, opts );
%}
fprintf('\nTraining BAT Network ....');
[Mynet_bat diagnostic]  = trainNetwork(...
    trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,3), Mynet, opts );

figure; plot(log10(diagnostic.TrainingLoss)); drawnow; title('loss')
%{    
fprintf('\nTraining CBV Network ....');
Mynet_cbv = trainNetwork(trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,2), Mynet, opts );

fprintf('\nTraining BAT Network ....');
Mynet_bat = trainNetwork(trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,3), Mynet, opts );

fprintf('\nTraining R1 Network ....');
Mynet_r1 = trainNetwork(trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,4), Mynet, opts );
%}
toc

%save NN_stuff Mynet aq_parms Nframes 

%%

fprintf('\nTesting Network with new synthetic data (linear)....');
% put the parms in a vector:
Ntests = 100

% test the network with different parameters
% with combinations from the training data
% recall ....
% trainingParms = [fvals cbvavals batvals  r1vals flipvals alpha_tivals];

idx = randperm(Ntrain); idx = idx(1:Ntests);
%tfvals = trainingParms(idx, 1);
tfvals = trainingParms(idx, 1) ;

idx = randperm(Ntrain); idx = idx(1:Ntests);
tcbvavals = trainingParms(idx, 2);

idx = randperm(Ntrain); idx = idx(1:Ntests);
tbatvals = trainingParms(idx, 3) ;

idx = randperm(Ntrain); idx = idx(1:Ntests);
tr1vals = trainingParms(idx, 4);

idx = randperm(Ntrain); idx = idx(1:Ntests);
tflipvals = trainingParms(idx,5);

mtis0vals = 1;

% put them together:
testparms = [tfvals tcbvavals  tbatvals  tr1vals tflipvals];
estimates = zeros(size(testparms));

test_parms = parms;  % initalize to the default values
n=1;
mse=[];
n=1;
%
test_noise = 0; %1e-3;
for n=1:Ntests
    
    truth = testparms(n,:) ;
    
    test_parms.f = truth(1);
    test_parms.cbva = truth(2);
    test_parms.bat = truth(3);
    test_parms.r1tis = truth(4);
    % test_parms.flip = truth(5);
    
    test_data = gen_signals_vs_201215(test_parms , aq_parms, doFigs, doSub, dt);
    %test_data = gen_signals_vs_201020(test_parms , aq_parms, doFigs, doSub, dt);
    Nframes = length(test_data);
    test_data = test_data/test_data(1);
    test_data = reshape(test_data, Nframes, 1,1);
    test_data = test_data + randn(size(test_data))*test_noise;

    % plot(test_data);drawnow
    % est_cbf = predict(Mynet_cbf, test_data);
    % est_cbv = predict(Mynet_cbv, test_data)
    est_bat = predict(Mynet_bat, test_data);
    % est_R1 = predict(Mynet_r1, test_data);

    
    % calc error
    %prediction = [est_cbf est_cbv est_bat est_R1];
    %mse(n) = mean(sqrt((truth(1:4)  - prediction).^2) ./ truth(1:4));
    
    % scale the prediciotn back  (see above scaling vector)
    %estimates(n,1) = est_cbf;
    %estimates(n,2) = est_cbv;
    estimates(n,3) = est_bat;
    %estimates(n,4) = est_R1;    
    
end

mytitles={
    'perfusion', 'CBva', 'BAT', 'R1', 'Flip'
    };

figure
n=3
%for n=1
%    subplot(3,2,n)
    plot(testparms(:,n), estimates(:,n), 'x')
    r = corrcoef(testparms(:,n), estimates(:,n));
    X = [testparms(:,n) ones(Ntests,1)];
    betas = pinv(X)*estimates(:,n);
    axis([min(testparms(:,n))  1.1*max(testparms(:,n)) min(testparms(:,n))  1.1*max(testparms(:,n))])
    title( (mytitles(n)))
    xlabel(sprintf('r = %0.2f, slope = %0.2f, int=%0.3f', r(1,2), betas(1), betas(2)))
%end

% subplot(3,2,6)
% plot([1:Ntests], mse)
% xlabel('% change in parameters')
% ylabel('MSRE')

return


