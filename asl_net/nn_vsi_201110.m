clear all
%
addpath  /home/hernan/matlab/flow/fingerprint
read_parms_file = 0
dt = 1e-4;
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
    aq_parms.labelcontrol           = ones(Nframes,1);
    aq_parms.labelcontrol(2:2:end)  = 0;
%    aq_parms.labelcontrol(1:6)        = -1;
    
    aq_parms.t_aq                   = 0.5 *ones(Nframes,1);
    aq_parms.ArtSup_delay           = 0.1 *ones(Nframes,1); % delay between AS pulse and acqusition
    
    aq_parms.doArtSup               = ones(Nframes,1);
%     aq_parms.doArtSup(5:5:end/2)    = 0;
%     aq_parms.doArtSup(1:4)          = -1;
    
    aq_parms.t_adjusts              = 1.5 *ones(Nframes,1);
    aq_parms.t_tag                  = zeros(Nframes,1);
    aq_parms.t_delays               =  0.20 + 1.5*abs(chirp(linspace(0,1,Nframes), 0, 1,3.5))';
    aq_parms.t_delays(2:2:end)      = aq_parms.t_delays(1:2:end);
   

end

% default parameters
parms.f         = 60/6000;
parms.Mtis0     = 1;
parms.cbva      = 0.02;
parms.bat       = 0.5;
parms.r1tis     = 1/1.4;
parms.flip      = deg2rad(90);
parms.alpha_ai  = 0.8; % arterial inversion efficiency
parms.alpha_ti  = 0.75; % tissue inversion efficiency
parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression



% sanity check: show the signal for one case
doSub = 0;
doFigs = 0;

subplot(211)
entry1 = abs(gen_signals_vs_201020(parms, aq_parms, doFigs,doSub, dt));
plot(entry1/entry1(1)); 
hold on

%parms.r1tis = 0.7*parms.r1tis;
parms.f = 30/6000;
entry2 = abs(gen_signals_vs_201020(parms, aq_parms, doFigs,doSub, dt));
plot(entry2/entry2(1)); 
legend('f = 60', 'f = 30')
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
% Ntrain = 1e6;
Nparms = 6;

mtis0vals = 1;
fvals =    single( rand(Ntrain, 1) * 0.02);     % ml/ml/s
cbvavals = single( rand(Ntrain, 1) * 0.1 );     % fraction
batvals = single( (rand(Ntrain, 1) * 0.5 + 0.05) );     % s
r1vals =   single( (rand(Ntrain, 1)* 3 + 0.3));     % 1/s

flipvals = single( pi/2*(0.995 + 0.01*(rand(Ntrain, 1))));     % rad : 90 +/- 10% 
%flipvals(:) = pi/2;

alpha_tivals  = single(rand(Ntrain, 1)*0.1 + 0.7 );  % tissue inversion efficiency  <--- reduced variability.
alpha_tivals(:)  = single( 0.75);  % tissue inversion efficiency  <--- reduced variability.


% allocate space and store training parms into single matrix
Nframes = length(entry1);

Ncombinations = Ntrain;

trainingData = single(zeros(Ncombinations, Nframes));

%trainingParms = single(zeros(Ncombinations, Nparms));

trainingParms = [fvals cbvavals batvals  r1vals flipvals alpha_tivals];

% keep as much memory clear as possible
clear fvals cbvavalsbatvals  r1vals flipvals

%  synth. training data
fprintf('\nGenerating training data ....');
tic
doFigs = 0;
doSub = 0;

parfor n=1:Ncombinations
    
    % fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    parms = struct();

    parms.mtis0     =   1;
    parms.alpha_ai  = 0.6; % arterial inversion efficiency
    parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression
    parms.f =       trainingParms(n, 1) ;
    parms.cbva =    trainingParms(n, 2) ;
    parms.bat =     trainingParms(n, 3) ;
    parms.r1tis =   trainingParms(n, 4) ;
    parms.flip =    trainingParms(n, 5) ;
    parms.alpha_ti= trainingParms(n, 6) ; % tissue inversion efficiency

    entry = abs(gen_signals_vs_201020(parms, aq_parms, doFigs, doSub, dt));
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
c1 = convolution2dLayer([round(Nframes/10) 1], 2,'stride',1);
r1 = reluLayer();
f0 = fullyConnectedLayer(Nframes);
bn = batchNormalizationLayer;

f1 = fullyConnectedLayer(Nframes);
f2 = fullyConnectedLayer(Nframes/2);
f3 = fullyConnectedLayer(round(Nframes/4));
f4 = fullyConnectedLayer(1);
r6 = regressionLayer;

%Mynet=[inputLayer; f1; r1; f2; r1; f3; r1; f4; r1;  r6]

% Try including batch normalization:
 Mynet=[inputLayer;   f1; r1; bn;   f2; r1; bn;   f3; r1; bn; f4;  r6]

delete(gcp)

opts = trainingOptions('sgdm', 'InitialLearnRate', 1e-4, ...
    'MaxEpochs', 4, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'parallel', ...
    'WorkerLoad', ones([40 1]))

save NN_stuff Mynet aq_parms Nframes  trainingData trainingParms

tic
fprintf('\nTraining CBF Network ....');
noise_leval = 0.005;
Mynet_cbf = trainNetwork(trainingData + randn(size(trainingData))*noise_level , ...
    trainingParms(:,1), Mynet, opts );

%{    
fprintf('\nTraining BAT Network ....');
Mynet_bat = trainNetwork(trainingData, trainingParms(:,3), Mynet, opts )

fprintf('\nTraining CBV Network ....');
Mynet_cbv = trainNetwork(trainingData, trainingParms(:,2), Mynet, opts )

fprintf('\nTraining R1 Network ....');
Mynet_r1 = trainNetwork(trainingData, trainingParms(:,4), Mynet, opts )
%}
toc

%save NN_stuff Mynet aq_parms Nframes 

%
fprintf('\nTesting Network with new synthetic data (linear)....');
% put the parms in a vector:
Ntests = 50

% test the network with different parameters
mtis0vals = 1;
fvals    = single( 0.02*randperm(Ntests)/Ntests )';            % ml/ml/s
cbvavals = single( 0.05*randperm(Ntests)/Ntests )';             % fraction
batvals  = single( 0.5 *randperm(Ntests)/Ntests )' ;     % s
r1vals   = single( 2.5*(randperm(Ntests)/ Ntests) )';       % 1/s
flipvals = single( (0.05*randn(size(r1vals)) +1)* pi/2 );     % rad ~50 to 100 degrees

testparms = [fvals cbvavals  batvals  r1vals flipvals];
estimates = zeros(size(testparms));

test_parms = struct();
n=1
mse=[];
n=1;
%
for n=1:Ntests
    
    truth = testparms(n,:) ;
    
    test_parms.f = truth(1);
    test_parms.cbva = truth(2);
    test_parms.bat = truth(3);
    test_parms.r1tis = truth(4);
    test_parms.flip = truth(5);
    test_parms.mtis0 = 1;
    test_parms.alpha_ai  = 0.8; % arterial inversion efficiency
    test_parms.alpha_ti  = 0.75; % tissue inversion efficiency
    test_parms.alpha_ts  = 0.17; % T2 weighting in tissue due to arterial suppression

    test_data = gen_signals_vs_201020(test_parms , aq_parms, doFigs, doSub, dt);
    Nframes = length(test_data);
    test_data = test_data/test_data(1);
    test_data = reshape(test_data, Nframes, 1,1);
    plot(test_data);drawnow
    est_cbf = predict(Mynet_cbf, test_data)
   % est_cbv = predict(Mynet_cbv, test_data)
   % est_bat = predict(Mynet_bat, test_data);
   % est_R1 = predict(Mynet_r1, test_data);

    
    % calc error
    %prediction = [est_cbf est_cbv est_bat est_R1];
    %mse(n) = mean(sqrt((truth(1:4)  - prediction).^2) ./ truth(1:4));
    
    % scale the prediciotn back  (see above scaling vector)
    estimates(n,1) = est_cbf;
    %estimates(n,2) = est_cbv;
    %estimates(n,3) = est_bat;
    %estimates(n,4) = est_R1;    
    
end

mytitles={
    'perfusion', 'CBva', 'BAT', 'R1', 'Flip'
    };

figure
for n=1:1
    subplot(3,2,n)
    plot(testparms(:,n), estimates(:,n), 'x')
    %axis([min(testparms(:,n))  max(testparms(:,n)) min(testparms(:,n))  max(testparms(:,n))])
    title( (mytitles(n)))
end

% subplot(3,2,6)
% plot([1:Ntests], mse)
% xlabel('% change in parameters')
% ylabel('MSRE')

return


