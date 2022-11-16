Npoints = 100;
Ntrain = 10000;
t=linspace(-1,1,Npoints);
trainingParms = 1*rand(Ntrain,2);
trainingData = zeros( Npoints, Ntrain);

for n=1:Ntrain
    a = trainingParms(n,1);
    b = trainingParms(n,2);
    trainingData(:,n) = (a + cos(b*t))';
end
    
trainingData = reshape(trainingData,Npoints, 1,1,Ntrain);

testData = 0.56 + cos(0.67*t);
testData = testData + 0.1*randn(size(testData));
testData = reshape(testData, Npoints, 1,1);

inputLayer = imageInputLayer([Npoints 1]);
c1 = convolution2dLayer([round(Npoints/10) 1], 2,'stride',1);
f1 = fullyConnectedLayer(50);
r1 = reluLayer();
f2 = fullyConnectedLayer(25);
f3 = fullyConnectedLayer(12);
f4 = fullyConnectedLayer(6);
f5 = fullyConnectedLayer(2);
r6 = regressionLayer;

Mynet=[inputLayer; f1; r1; f2; r1; f5; r1; r6]

opts = trainingOptions('sgdm', 'InitialLearnRate', 1e-3, 'MaxEpochs', 500)

Mynet = trainNetwork(trainingData, trainingParms, Mynet, opts )

prediction = predict(Mynet, testData)


%



