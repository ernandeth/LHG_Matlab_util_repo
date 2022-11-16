fprintf('\nloading Network and timing parms....');
load timing_parms
load(netName);

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
%bat1vals = single( linspace(2.7, 0.2,  Ntests))' ;     % s

r1vals =   single( linspace(0.3, 2.8, Ntests) )';       % 1/s
flipvals = single( linspace(0.87, 1.74,Ntests) )';     % rad ~50 to 100 degrees

testparms = [fvals cbvavals kforvals bat1vals  r1vals flipvals];
estimates = zeros(size(testparms));

test_parms = struct();
n=1
mse_lin=[];
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
    
    % adding noise to the testing data
    test_data = test_data + 0.001*randn(size(test_data));
    %test_data = test_data / test_data(1);   % Normalize to the first data point. 
    test_data = test_data/mean(test_data);

    test_data = reshape(test_data, Npoints, 1,1);
    prediction = predict(Mynet, test_data);
    
    
    % calc error
    truth = truth .* scale_vector;
    mse_lin(n) = mean(sqrt((truth  - prediction).^2) ./ truth);
    
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
rhos_lin = zeros(1,6);
for n=1:6
    subplot(3,2,n)
    plot(testparms(:,n), estimates(:,n), '.'); hold on;
    plot([0 1], [0 1], 'r'); hold off
    tmp = corrcoef(testparms(:,n), estimates(:,n));
    rhos_lin(n) = tmp(1,2);
    ylabel(sprintf('RHO=%f', rhos_lin(n)))
    
    axis([min(testparms(:,n))  max(testparms(:,n)) min(testparms(:,n)) max(testparms(:,n))])
    title( (mytitles(n)))
    axis square
end
set(gcf,'Name', netName)

%{
figure
plot([1:Ntests], mse_lin)
xlabel('% change in parameters')
ylabel('MSRE')
%}

%% a random combination test
fprintf('\n Testing network with random Parameters ....');
dofigs=0;
doSub=0;

Ntests = 5000;

mtis0vals = 1;
Dispvals = [20];
fvals =    single( rand(Ntests, 1) * 0.02);            % ml/ml/s
cbvavals = single( rand(Ntests, 1) * 0.1);             % fraction
kforvals = single( rand(Ntests, 1) * 0.5);             % 1/s
bat1vals = single( (rand(Ntests, 1) * 2.5 + 0.2));     % s
%bat2vals = single( (rand(Ntests, 1) * 2.5 + 0.2));     % s
r1vals =   single( (rand(Ntests, 1)*2.5 + 0.3));       % 1/s
flipvals = single( (rand(Ntests, 1)*0.87 + 0.87));     % rad ~50 to 100 degrees

mse_rand = [];
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
    
    % adding noise to the testing data
    test_data = test_data + 0.001*randn(size(test_data));
    
    % scale to the first data point AFTER noise addition
    %test_data = test_data/test_data(1);
    test_data = test_data/mean(test_data);
    
    
    test_data = reshape(test_data, Npoints, 1,1);
    prediction = predict(Mynet, test_data);
    
    % But in physical units, we need to scale the prediciotn vector  (see above scaling vector)
    prediction= prediction ./ scale_vector;
    estimates(n,:) = prediction;
    
    % calc error (MSE)
    mse_rand(n) = mean(sqrt((truth  - prediction).^2) ./ truth);
    
end

mytitles={
    'perfusion', 'CBva', 'kfor', 'BAT', 'R1', 'Flip'
    }
figure

rhos_rand = zeros(1,6);
for n=1:6
    subplot(3,2,n)
    plot(testparms(:,n), estimates(:,n), '.'); hold on;
    plot([0 1], [0 1], 'r'); hold off
    tmp = corrcoef(testparms(:,n), estimates(:,n));
    rhos_rand(n) = tmp(1,2);
    ylabel(sprintf('RHO=%f', rhos_rand(n)))
    
    axis([min(testparms(:,n))  max(testparms(:,n)) min(testparms(:,n)) max(testparms(:,n))])
    title( (mytitles(n)))
    axis square
end
set(gcf,'Name', netName)

% testing the linearity of the perfusion fit

y = estimates(:,1);
x = [ones(Ntests,1) testparms(:,1) testparms(:,1).^2];
betas = pinv(x)*y
subplot(3,2,1)
hold on
plot(x(:,2), x*betas,'g-')
hold off
%{
figure
plot(mse_rand,'x')
xlabel('% change in parameters')
ylabel('MSRE')
%}
