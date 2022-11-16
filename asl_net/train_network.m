%
%clear all
%
addpath  /home/hernan/matlab/flow/fingerprint
netSize = 9;
makeNewTrainingData = 0;

%% create ASL sequence timing
Npoints = 300;

% scaling the output of the network so that the
% norm of the vector containing the output (estimates) gives the right
% amount of weight to each parameter when calculating the error.
% f(~0.01), CBV(~0.01), kfor(~0.5), BAT(~1.2)  R1(~1.0)  flip(~1)
% scale_vector = [800 100 1 5  1 1 ]
% scale_vector = [300 100 1 5  1 1 ]
% scale_vector = [600 100 1 5  1 1 ]
% scale_vector = [400 100 1 5  1 1 ]
scale_vector = [500 100 1 5  1 1 ]

typical = [0.01 0.01 0.3 1.2 1 1]
scale_vector .* typical


if makeNewTrainingData
    timing_parms.t_tag = 2*abs(sinc(linspace(2, 0, Npoints))) + 0.05 ;
    timing_parms.t_delay =  0.05 * ones(1,Npoints) ;
    
    %timing_parms.t_tag = 2*abs( (linspace(0.1,1, Npoints)) .* cos(linspace(0,3*pi,Npoints))) + 0.05 ;
    %timing_parms.t_tag = 2*abs(sinc(linspace(-3, 3, Npoints))) + 0.05 ;
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
    load minLinInterp360LH350_TR_fix.mat
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
        'kfor', 0.02, ...
        'r1tis', 0.8,...
        'flip', pi/4, ...
        'Disp', 20);
    
    
    % show the signal for one case
    doSub = 0;
    dofigs = 1;
    
    entry = gen_signals_180320(parms0, timing_parms, dofigs, doSub);
    drawnow
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
    % fvals =    single( rand(Ntrain, 1) * 0.02       * scale_vector(1));     % ml/ml/s
    % fvals =    single( rand(Ntrain, 1) * 0.03       * scale_vector(1));     % ml/ml/s
    % fvals =    single( rand(Ntrain, 1) * 0.04       * scale_vector(1));     % ml/ml/s
    % fvals = single ([abs(0.25*randn(Ntrain/2,1)); 1-abs(0.25*randn(Ntrain/2,1))]*0.03 * scale_vector(1) );
    % fvals = single ([abs(0.25*randn(Ntrain/2,1)); 1-abs(0.25*randn(Ntrain/2,1))]*0.025 * scale_vector(1) );
    fvals = single (abs([abs(0.35*randn(Ntrain/2,1)); 1-abs(0.45*randn(Ntrain/2,1))])*0.03);
    fvals(fvals>0.03)=0.03;
    fvals = fvals * scale_vector(1);
    
    cbvavals = single( rand(Ntrain, 1) * 0.1        * scale_vector(2));     % fraction
    kforvals = single( rand(Ntrain, 1) * 0.5        * scale_vector(3));     % 1/s
    bat1vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(4));     % s
    %bat2vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(5));     % s
    r1vals =   single( (rand(Ntrain, 1)*2.5 + 0.3)  * scale_vector(5));     % 1/s
    flipvals = single( (rand(Ntrain, 1)*0.87 + 0.87)* scale_vector(6));     % rad ~50 to 100 degrees
    %
    %     mtis0vals = 1;
    %     Dispvals = [20];
    %     fvals =    single( rand(Ntrain, 1) * 0.025       * scale_vector(1));     % ml/ml/s
    %     cbvavals = single( rand(Ntrain, 1) * 0.1        * scale_vector(2));     % fraction
    %     kforvals = single( rand(Ntrain, 1) * 0.5        * scale_vector(3));     % 1/s
    %     bat1vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(4));     % s
    %     %bat2vals = single( (rand(Ntrain, 1)*2.5 + 0.2)  * scale_vector(5));     % s
    %     r1vals =   single( (rand(Ntrain, 1)*2.5 + 0.3)  * scale_vector(5));     % 1/s
    %     flipvals = single( (rand(Ntrain, 1)*pi/2 + 0.20)* scale_vector(6));     % rad ~50 to 100 degrees
    %
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
    % adding noise to the TRAINING data
    trainingData = trainingData + 0.005*randn(size(trainingData));
    
    % scaling the training data AFTER the noise
    for n=1:Ncombinations
        %trainingData(n,:) = trainingData(n,:) / trainingData(n,1);
        trainingData(n,:) = trainingData(n,:) / mean(trainingData(n,:));
    end
    
    % reshape it so that it fits in the Matlab NN framework
    trainingData = reshape(trainingData', Npoints, 1,1,Ncombinations);
    %
    whos
    save -v7.3 trainingData.mat trainingData
    save -v7.3 trainingParms.mat trainingParms
    save timing_parms.mat timing_parms
    
    % files for the scanner:
    !mkdir new_schedule
    tmp = timing_parms.t_tag;
    save new_schedule/t_tags.txt tmp -ascii
    tmp = timing_parms.t_adjust;
    save new_schedule/t_adjusts.txt tmp -ascii
    tmp = timing_parms.t_delay;
    save new_schedule/t_delays.txt tmp -ascii
    tmp = timing_parms.isLabel;
    save new_schedule/labelcontrol.txt tmp -ascii
    
    %
    % try differencing the data:
    % D_trainingData = trainingData(2:2:end, :,:,:) - trainingData( 1:2:end, :,:,:);
    % Npoints = Npoints/2;
    %
else
    
    if ~exist('trainingData')
        fprintf('\nLoading training data from files ....');
        
        load trainingData.mat
        load trainingParms.mat
        load timing_parms.mat
    end
    
    % noise test... just adding more noise to the same training data!
    % trainingData = trainingData + 0.015*randn(size(trainingData));
    
    Npoints = size(trainingData,1);
    
    
    
end
%%
%}
fprintf('\nGenerating Network ....');

% a selection of layers for the network
inputLayer = imageInputLayer([Npoints 1]);
seqIn = sequenceInputLayer(1);
c1 = convolution2dLayer([round(Npoints/20) 1], 2,'stride',1);
r1 = reluLayer();
f150 = fullyConnectedLayer(150);
f100 = fullyConnectedLayer(100);
f050 = fullyConnectedLayer(50);
f025 = fullyConnectedLayer(25);
f012 = fullyConnectedLayer(12);
f006 = fullyConnectedLayer(6);
r006 = regressionLayer;


clear Mynet

% configure the network here with the above layers
switch(netSize)
    case 1
        Mynet= [inputLayer; f150; r1; f150; r1; f006; r1; r006]
    case 2
        Mynet= [inputLayer; f100; r1; f100; r1; f006; r1; r006]
    case 3  % this one seems the best
        Mynet= [inputLayer; f050; r1; f050; r1; f006; r1; r006]
    case 4
        Mynet= [inputLayer; f100; r1; f050; r1; f025; r1; f006; r1; r006]
    case 5
        Mynet= [inputLayer; f050; r1; f050; r1; f050; r1; f006; r1; r006]
    case 6
        Mynet= [inputLayer; c1;   r1; f100; r1; f100; r1; f006; r1; r006]
    case 7
        Mynet= [inputLayer; f050; r1; f025; r1; f012; r1; f006; r1; r006]
    case 8
        Mynet= [inputLayer; f050; r1; f050; r1; f012; r1; f006; r1; r006]
    case 9
        Mynet= [inputLayer; f050; r1; f025; r1;  f006; r1; r006]
        
end



delete(gcp)

opts = trainingOptions('sgdm', 'InitialLearnRate', 1e-3, ...
    'MaxEpochs', 3, ...
    'Plots', 'none', ...
    'ExecutionEnvironment', 'parallel', ...
    'WorkerLoad', ones([40 1]))


tic
fprintf('\nTraining Network ....');
Mynet = trainNetwork(trainingData, trainingParms, Mynet, opts )
toc

switch(netSize)
    case 1
        save NN_2layer_150_150.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_2layer_150_150.mat'
    case 2
        save NN_2layer_100_100.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_2layer_100_100.mat'
        
    case 3
        save NN_2layer_050_050.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_2layer_050_050.mat'
        
    case 4
        save NN_3layer_100_050_025.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_3layer_100_050_025.mat'
        
    case 5
        save NN_3layer_050_050_50.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_3layer_050_050_50.mat'
        
    case 6
        save NN_2layer_conv_100_100.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_2layer_conv_100_100.mat'
        
    case 7
        save NN_3layer_050_025_012.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_3layer_050_025_012.mat'
        
    case 8
        save NN_3layer_050_050_012.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_3layer_050_050_012.mat'
        
    case 9
        save NN_2layer_050_025.mat Mynet timing_parms Npoints scale_vector
        netName = 'NN_2layer_050_025.mat'
        
end

return
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
%}