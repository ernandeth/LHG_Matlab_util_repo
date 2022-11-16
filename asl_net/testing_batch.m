netNames = ...
{
'NN_2layer_150_150.mat',
'NN_2layer_100_100.mat',
'NN_2layer_050_050.mat',
'NN_3layer_100_050_025.mat',
'NN_3layer_050_050_050.mat',
'NN_2layer_conv_100_100.mat', 
'NN_3layer_050_025_012.mat',
'NN_3layer_050_050_012.mat',
'NN_2layer_050_025.mat'};

close all

all_mse_lin=[];
all_mse_rand = [];
all_rhos_rand = [];
all_rhos_lin = [];

plotnum=1;
for nn=1:length(netNames)
   netName = netNames{nn}
   
   test_network;
    
    all_mse_lin = [all_mse_lin; mse_lin]; 
    all_mse_rand = [all_mse_rand; mse_rand]; 

    all_rhos_lin = [all_rhos_lin; rhos_lin];
    all_rhos_rand = [all_rhos_rand; rhos_rand];
    
    figure(300)
    subplot(3,3,plotnum)
    plot(mse_rand)
    title(netName, 'Interpreter', 'none')
    axis([1 500 0 5])

    figure(301)
    subplot(3,3,plotnum)
    plot(mse_lin)
    title(netName, 'Interpreter', 'none')
    axis([1 50 0 0.5])

    plotnum = plotnum+1;
end


figure(300)
set(gcf,  'Name', 'individual RMSE (randomized)')
figure(301)
set(gcf,  'Name', 'individual RMSE (linearized)')

counter = all_mse_rand;
threshold = 1
counter(counter <= threshold) = 1;
counter(counter > threshold) = 0;
netNames

eval_table = [all_rhos_lin all_rhos_rand sum(all_mse_rand,2),    sum(counter,2)]
