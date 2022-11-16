clear all
%
addpath  /home/hernan/matlab/flow/fingerprint

%% create ASL sequence timing
Npoints = 300;

% scaling the output of the network so that the
% norm of the vector containing the output (estimates) gives the right
% amount of weight to each parameter when calculating the error.
% f(~0.01), CBV(~0.01), kfor(~0.5), BAT(~1.2)  R1(~1.0)  flip(~1)
scale_vector = [500 100 1 5  1 1 ]
typical = [0.01 0.01 0.3 1.2 1 1]
scale_vector .* typical

timing_parms.t_tag = 2*abs(sinc(linspace(-2,0, Npoints))) + 0.05 ;
timing_parms.t_delay =  0.05 * ones(1,Npoints) ;

%timing_parms.t_tag = 2*abs( (linspace(1, 0.1,Npoints)) .* cos(linspace(0,3*pi,Npoints))) + 0.05 ;
%timing_parms.t_tag = 0.5 + 1.5*rand(1,Npoints)  ;
%timing_parms.t_delay = 0.05 + 0.4*abs(sin(linspace(0, 1*pi,Npoints))) ;

timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;

timing_parms.isLabel = round(rand(1,Npoints))  ;
supercontrols = randperm(Npoints);
supercontrols=supercontrols(1:round(Npoints/4));
timing_parms.isLabel(supercontrols) = -1;

timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;

% testing out anish's schedule:
%{
load minPoly360LH500_20_TR_fix.mat
timing_parms = timing_parms_min;
Npoints = length(timing_parms.t_tag);
%}

%
%synth. training parameters (targets)
parms0 = struct( ...
    'mtis0', 1,...
    'f', 0.01 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0.5, ...
    'r1tis', 0.8,...
    'flip', pi/4, ...
    'Disp', 20);


% show the signal for one case
doSub = 0;
dofigs = 1;

entry = gen_signals_180320(parms0, timing_parms, dofigs, doSub);
dofigs = 0;

%%

fprintf('\nGenerating training Parameters ....');
% scale the training parameters so that the weights can be adjusted in the
% RMSE calculation at the end of the network  during backpropagation
n = 1;
Ntrain =  1.2e7;      % memory used for training parms: Ntrain x 7 x 4 Bytes.
% memory for training data is : Ntrain x Npoints x 4 Bytes
% This means if Ntrain = 10^7
% training parms size is 280 MB
% training data size is 12,000 MB

% again, for Anish's test: 
% Ntrain = 1e6;

mtis0vals = 1;
Dispvals = [20];
fvals =    single( rand(Ntrain, 1) * 0.02       * scale_vector(1));     % ml/ml/s
cbvavals = single( rand(Ntrain, 1) * 0.1        * scale_vector(2));     % fraction
kforvals = single( rand(Ntrain, 1) * 0.5        * scale_vector(3));     % 1/s
bat1vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(4));     % s
%bat2vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(5));     % s
r1vals =   single( (rand(Ntrain, 1)*2.5 + 0.3)  * scale_vector(5));     % 1/s
flipvals = single( (rand(Ntrain, 1)*0.87 + 0.87)* scale_vector(6));     % rad ~50 to 100 degrees



% allocate space and store training parms into single matrix
Npoints = length(entry);

Ncombinations = Ntrain;

trainingData = single(zeros(Ncombinations, Npoints));

trainingParms = single(zeros(Ncombinations, 6));
trainingParms = [fvals cbvavals kforvals bat1vals  r1vals flipvals];

% keep as much memory clear as possible
clear fvals cbvavals kforvals bat1vals  r1vals flipvals

%  synth. training data
fprintf('\nGenerating training data ....');

parfor n=1:Ncombinations
    
    % fprintf('\rGenerating entry .... %d    of  %d   ',n, Ncombinations);
    parms = struct();
    
    % undo the scaling in order to generate the correct signals
    parms.f =       trainingParms(n, 1) ./ scale_vector(1);
    parms.cbva =    trainingParms(n, 2) ./ scale_vector(2);
    parms.kfor =    trainingParms(n, 3) ./ scale_vector(3)
    parms.bat =     trainingParms(n, 4) ./ scale_vector(4);
    parms.bat2 =    0; % gets ignored
    parms.r1tis =   trainingParms(n, 5) ./ scale_vector(5);
    parms.flip =    trainingParms(n, 6) ./ scale_vector(6);
    parms.mtis0=    1;
    parms. Disp =   20;
    
    entry = gen_signals_180320(parms, timing_parms, dofigs, doSub);
    
    trainingData(n,:) = entry;
    
end

trainingData = reshape(trainingData', Npoints, 1,1,Ncombinations);
whos

%%
% try differencing the data:
% D_trainingData = trainingData(2:2:end, :,:,:) - trainingData( 1:2:end, :,:,:);
% Npoints = Npoints/2;

%
fprintf('\nGenerating Network ....');

inputLayer = imageInputLayer([Npoints 1]);
c1 = convolution2dLayer([round(Npoints/10) 1], 2,'stride',1);
r1 = reluLayer();
f0 = fullyConnectedLayer(100);
f1 = fullyConnectedLayer(50);
f2 = fullyConnectedLayer(25);
f3 = fullyConnectedLayer(12);
f4 = fullyConnectedLayer(6);
r6 = regressionLayer;

Mynet=[inputLayer; f1; r1; f2; r1; f4; r1;  r6]

delete(gcp)

opts = trainingOptions('sgdm', 'InitialLearnRate', 1e-3, ...
    'MaxEpochs', 4, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'parallel', ...
    'WorkerLoad', ones([40 1]))

save NN_stuff Mynet timing_parms Npoints scale_vector trainingData trainingParms

tic
fprintf('\nTraining Network ....');
Mynet = trainNetwork(trainingData, trainingParms, Mynet, opts )
toc

save NN_stuff Mynet timing_parms Npoints scale_vector

%%
fprintf('\nTesting Network with new synthetic data (linear)....');
% put the parms in a vector:
Ntests = 50

% test the network with different parameters
mtis0vals = 1;
Dispvals = [20];
fvals =    single( linspace(0, 0.02,Ntests) )';            % ml/ml/s
cbvavals = single( linspace(0, 0.1, Ntests) )';             % fraction
kforvals = single( linspace(0, 0.5, Ntests) )';             % 1/s
bat1vals = single( linspace(0.2, 2.7, Ntests))' ;     % s
r1vals =   single( linspace(0.3, 2.8, Ntests) )';       % 1/s
flipvals = single( linspace(0.87, 1.74,Ntests) )';     % rad ~50 to 100 degrees

testparms = [fvals cbvavals kforvals bat1vals  r1vals flipvals];
estimates = zeros(size(testparms));

test_parms = struct();
n=1
mse=[];
n=1;

for n=1:Ntests
    
    truth = testparms(n,:) ;
    
    test_parms.f = truth(1);
    test_parms.cbva = truth(2);
    test_parms.kfor = truth(3);
    test_parms.bat = truth(4);
    test_parms.bat2 = 0;
    test_parms.r1tis = truth(5);
    test_parms.flip = truth(6);
    test_parms.mtis0 = 1;
    test_parms.Disp = 20;
    
    test_data = gen_signals_180320(test_parms , timing_parms, 0, 0);
    test_data = reshape(test_data, Npoints, 1,1);
    prediction = predict(Mynet, test_data);
    
    
    % calc error
    truth = truth .* scale_vector;
    mse(n) = mean(sqrt((truth  - prediction).^2) ./ truth);
    
    % scale the prediciotn back  (see above scaling vector)
    prediction= prediction ./ scale_vector;
    estimates(n,:) = prediction;
end

test_data = gen_signals_180320(test_parms , timing_parms, 1, 1);
%hold on, plot( test_data)



mytitles={
    'perfusion', 'CBva', 'kfor', 'BAT', 'R1', 'Flip'
    }
figure
for n=1:6
    subplot(4,2,n)
    plot(testparms(:,n), estimates(:,n), 'x')
    axis([min(testparms(:,n))  max(testparms(:,n)) min(testparms(:,n)) max(testparms(:,n))])
    title( (mytitles(n)))
end

subplot(4,2,8)
plot([1:Ntests], mse)
xlabel('% change in parameters')
ylabel('MSRE')


%% a random combination test
fprintf('\n Testing network with random Parameters ....');

Ntests = 1000;

mtis0vals = 1;
Dispvals = [20];
fvals =    single( rand(Ntests, 1) * 0.02);            % ml/ml/s
cbvavals = single( rand(Ntests, 1) * 0.1);             % fraction
kforvals = single( rand(Ntests, 1) * 0.5);             % 1/s
bat1vals = single( (rand(Ntests, 1) * 2.5 + 0.2));     % s
%bat2vals = single( (rand(Ntests, 1) * 2.5 + 0.2));     % s
r1vals =   single( (rand(Ntests, 1)*2.5 + 0.3));       % 1/s
flipvals = single( (rand(Ntests, 1)*0.87 + 0.87));     % rad ~50 to 100 degrees


testparms = [fvals cbvavals kforvals bat1vals  r1vals flipvals];
estimates = zeros(size(testparms));
n=1;
for n=1:Ntests
    
    truth = testparms(n,:) ;
    
    test_parms.f = truth(1);
    test_parms.cbva = truth(2);
    test_parms.kfor = truth(3);
    test_parms.bat = truth(4);
    test_parms.bat2 = 0;
    test_parms.r1tis = truth(5);
    test_parms.flip = truth(6);
    test_parms.mtis0 = 1;
    test_parms.Disp = 20;
    
    test_data = gen_signals_180320(test_parms , timing_parms, dofigs, doSub);
    test_data = reshape(test_data, Npoints, 1,1);
    prediction = predict(Mynet, test_data);
    
    % But in physical units, we need to scale the prediciotn vector  (see above scaling vector)
    prediction= prediction ./ scale_vector;
    estimates(n,:) = prediction;
    
    % calc error (MSE)
    mse(n) = mean(sqrt((truth  - prediction).^2) ./ truth);
    
end

mytitles={
    'perfusion', 'CBva', 'kfor', 'BAT', 'R1', 'Flip'
    }
figure

for n=1:6
    subplot(4,2,n)
    plot(testparms(:,n), estimates(:,n), 'x')
    axis([min(testparms(:,n))  max(testparms(:,n)) min(testparms(:,n)) max(testparms(:,n))])
    title( (mytitles(n)))
end
subplot(4,2,8)
plot(mse,'x')
xlabel('% change in parameters')
ylabel('MSRE')

return




%% A quick test for gen_signal
%{

Npoints = 30;
durs = linspace(0.1, 3, Npoints);
timing_parms.t_tag =  durs ;
timing_parms.t_delay = 0.05 * ones(1,Npoints) ;
timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.isLabel =  ones(1,Npoints) ;
timing_parms.isLabel(2:2:end) =  0 ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;

parms0 = struct( ...
    'mtis0', 1,...
    'f', 0.01 , ...
    'cbva' ,  0.01 , ...
    'bat', 1.0,...
    'bat2', 1.2,...
    'kfor', 0.5, ...
    'r1tis', 0.8,...
    'flip', pi/4, ...
    'Disp', 20);
figure
sig=gen_signals_180320(parms0 , timing_parms, 1, 0);
figure
subplot(211)
plot(sig)
subplot(212)
plot(sig(1:2:end) - sig(2:2:end))

parms0.bat = 2.5;
sig=gen_signals_180320(parms0 , timing_parms, 0, 0);
subplot(211)
hold on
plot(sig)
hold off
subplot(212)
hold on
plot(sig(1:2:end) - sig(2:2:end))
hold off

%}
%%

% a test to see how much BAT affects the signal
%{
timing_parms.isLabel = round(rand(1,Npoints))  ;
supercontrols = randperm(Npoints);
supercontrols = supercontrols(1:round(Npoints/6));
timing_parms.isLabel(supercontrols) = -1;

timing_parms.t_delay = 0.05 + 0.4*abs(sin(linspace(0, 1*pi,Npoints))) ;
timing_parms.t_tag = 0.5 + 1.5*rand(1,Npoints)  ;

timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;
%}

test_parms = struct( ...
    'mtis0', 1,...
    'f', 0.001 , ...
    'cbva' ,  0.01 , ...
    'bat', 0.5,...
    'bat2', 1.2,...
    'kfor', 0.5, ...
    'r1tis', 0.8,...
    'flip', pi/4, ...
    'Disp', 150);

ref = gen_signals_180320(test_parms , timing_parms, 0, 0);
sd = zeros(1,10);
close
for n=1:10
    test_parms.bat = test_parms.bat + 0.2;
    test_data = gen_signals_180320(test_parms , timing_parms, 0, 0);
    
    subplot(211)
    plot(test_data); hold on
    subplot(212)
    plot(test_data-ref);hold on
    
    % what is the standard deviation of the difference introduced in teh
    % signal by change in BAT
    sd(n)=std(test_data-ref);
    %drawnow
    
end
% averaging over the range of BAT's that we tested:
mm = mean(abs(sd))

hold off
% figure
% plot(mm)
