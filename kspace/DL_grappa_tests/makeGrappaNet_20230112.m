function [myNet  rho] = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
% function [myNet test_correlation] = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
%
% (C) Luis Hernandez-Garcia @ University of Michigan
% hernan@umich.edu
%
% Learn a function that will map kspace  locations
% into signals for an individual data set
%
% based on work presented by Daniel Rueckert at Sedona 2023
%
% INPUTS:
% kx, ky, kz is the trajectory in k-space (cm*-1)
% data is the signal : (Npts_echo*Nechoes) x Ncoils
% dim3 is num. pixels in each dimension of IMAGE space
% CalRadius is the fraction of k space used for training region
%
% OUTPUT:
% returns a network that maps (kx,ky,kxz,coil) position to the corresponding signal.
%
%

% Remove Largest and smallest echos to try to reduced T2 effects
%{
len = size(data,1);
data = data(1:3*end/4, :);
kx = kx(1 : 3*end/4, :);
ky = ky(1 : 3*end/4, :);
kz = kz(1 : 3*end/4, :);
%}

% STEP 1: Set up a calibration region
% Determine a center region for calibration, assuming that this
% center region is well sampled.
% the calibration region is determined by CalRadius relative to the Rmax
R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
CalRadius = CalRadius * Kmax;    % convert from fraction to k space units.
% size of k-space voxels: (1/FOV) (distance between neighbors in the iso-tropic cartesian grid)
dkx = 2*Kmax/dim3(1);

% scale the kspace locations between -1 and 1
% Limit the size of the learning region 
% this is a way to avoid trouble regions
CalRegion = find( R>(Kmax*0.01) & R<CalRadius );
kx = kx(CalRegion) / Kmax;
ky = ky(CalRegion) / Kmax;
kz = kz(CalRegion) / Kmax;
data = data(CalRegion,:);

% reorganize the input data into a Nx4 matrix
Ncoils = size(data,2);
slen = length(data);

%
Inputs = zeros(slen*Ncoils, 4);  
for n=1:Ncoils
    data(:,n) = smooth(data(:,n),3);
    Inputs( [1:slen] + slen*(n-1), :) = [kx ky kz n*ones(slen,1)];
end

Outputs = data(:);

% debugging: see if we can learn a Gaussian
% Outputs = Inputs(:,4) .* exp(-(Inputs(:,1).^2 + Inputs(:,2).^2 + Inputs(:,3).^2 )/ (Kmax/20)); 

%{
%<--- let's try to do just one coil
Inputs = [ky ky kz];
Outputs = (data(:,10));  
%}


% split the data into training and testing 
test_inds = randperm(length(Outputs), ceil(length(Outputs)/20))';
train_inds = ones(size(Outputs));
train_inds(test_inds) = 0;
train_inds = shuffle(find(train_inds));
%train_inds = train_inds(1:end/2);

% convert into logs to even out the dynamic range
Outputs = log(abs(Outputs) + eps) .* exp(1i*angle(Outputs));
scale = max(abs(Outputs(:)));
Outputs = Outputs/scale;
% split into real and imaginary, add an offset and take the log
Outputs = [real(Outputs(:)) imag(Outputs(:))];
%{
Outputs(Outputs < 0) = -log(1 + abs(Outputs(Outputs < 0)));
Outputs(Outputs > 0) = log(1 + abs(Outputs(Outputs > 0)));
Outputs(Outputs == 0) = 0;
%}
%Outputs = Outputs.^0.2;

% weigh the data by their location: further data gets magnified
% we are trying to make the distribution a bit more even.
%{
W = sqrt(Inputs(:,1).^2 + Inputs(:,2).^2+ Inputs(:,3).^2);
Outputs = Outputs .* W .* W;
offset2 = mean(Outputs(:));
Outputs = Outputs - offset2;
%}

TrainDataIn = Inputs(train_inds,:);
TrainDataOut = Outputs(train_inds,:);

TestDataIn = Inputs(test_inds,:);
TestDataOut = Outputs(test_inds,:);

justafew=3:100:size(TestDataIn,1);

vData ={ TestDataIn(justafew,:), TestDataOut(justafew,:)};

hist(TrainDataOut(:),100);drawnow;

%
% CREATE NETWORK
%
% Goal: Learn a function that will map relative neighbor locations into
% signals at that k-space location
%

myNet = [...

featureInputLayer(4);

fullyConnectedLayer(4);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(8);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(8);
batchNormalizationLayer;
reluLayer;


fullyConnectedLayer(8);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(8);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(4);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(2);
regressionLayer
]

%{

layers = [...

featureInputLayer(4);

fullyConnectedLayer(10);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(10);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(2);
]

dlnet = dlnetwork(layers)
%}
% TRAIN NETWORK
%
fprintf('\nTraining Networks ...\n ');
options = trainingOptions('adam',...
    'Shuffle','every-epoch', ...
    'WorkerLoad', ones([40 1]), ...
    'L2Regularization',  1e-3, ... 
    'LearnRateSchedule', 'piecewise', ...  
    'MaxEpochs',5,...
    'Plots', 'training-progress',...
    LearnRateDropFactor = 0.5, LearnRateDropPeriod = 1);
    %'InitialLearnRate', 1e-4, ... %1e-4
    %'OutputNetwork', 'best-validation-loss', ...
    %'ValidationData', vData, ...
    %'LearnRateSchedule', 'piecewise', ...
    %'ValidationPatience',1e3, ...,
    %'MiniBatchSize', 50, ...  % reduce this (20) to avoid running out of memory (11.18.22)
    

% Now train this network once for each coil


[myNet, info] = trainNetwork(...
    TrainDataIn, ...
    TrainDataOut, ...
    myNet, ...
    options);


%fprintf('\nSaving networks ...')
%save GRAPPAnet.mat  myNet


%  TEST NETWORK
fprintf('\nTesting the network with data from calibration region...')

gpuDevice(1)
justafew=1:10:size(TestDataIn,1);

TestInput = TestDataIn(justafew,:);
Truth = TestDataOut(justafew,:);


% debugging - can we get it right, if we test the training data?
TestInput = TrainDataIn(justafew,:);
Truth = TrainDataOut(justafew,:);

%%%


est = zeros(size(Truth));


for n =1:length(justafew)
    tmp = predict(myNet , TestInput(n,:));
    est(n,:) = tmp;
end

inds = find(est==0);
rho = abs(corrcoef(Truth, est))

%
figure
plot(Truth(:,1), est(:,1),'.')
hold on
plot(Truth(:,2), est(:,2),'.r')

xlabel('Truth')
ylabel('Estimate')
legend(sprintf('Correlation %f', rho(2,1)))
line([-max(abs(Truth)) max(abs(Truth))], [ -max(abs(Truth))  max(abs(Truth))  ])
drawnow

return

