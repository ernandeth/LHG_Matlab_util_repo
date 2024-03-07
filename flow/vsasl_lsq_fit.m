function est = vsasl_lsq_fit(timeseries, aqparms, parms, showFit)

if nargin==0
    % Default values for testing:
    aqparms.del1 = 0.1;
    aqparms.del2= 0.3;
    aqparms.del3 = 0.05;
    aqparms.ro_times = [1:30]*0.05;

    parms.M0 = 1;
    parms.alpha = 0.86;
    parms.R1a = 1/1.6;  % sec-1
    parms.R1t = 1/1.4;  % sec-1
    parms.lambda = 0.9;
    showFit = 1;
end

f = 0.01;
CBVa = 0.015;
D = 0.2;

guess0 = [f , CBVa, D];
lb = [0 0 1e-10];
ub = [0.02, 0.05, 1];

if nargin==0
    % A test with synth data
    aqparms.del2= 0.3;
    timeseries1 = vsasl_fun(guess0, aqparms, parms, showFit, 0 );
    aqparms.del2= 1.5;
    timeseries2 = vsasl_fun(guess0, aqparms, parms, showFit, 0 );
clf
    subplot(211)
    plot(timeseries1); hold on
    plot(timeseries2);
    subplot(212)
    plot(timeseries2-timeseries1)
    
    timeseries = timeseries2 + 2e-4*randn(size(timeseries1));
    guess0 = guess0*5
    drawnow

end


opts = optimoptions(@lsqnonlin);
%opts.MaxIterations = 1000;
opts.SubproblemAlgorithm="cg";
% opts.FiniteDifferenceStepSize = 1e-8;
% %opts.MaxFunctionEvaluations = 1e3;
% opts.OptimalityTolerance = 1e-8;
% opts.FunctionTolerance = 1e-8;
% opts.StepTolerance = 1e-8;
opts.FiniteDifferenceType ='central';
%opts.Display = 'iter';
est=[];
%{
[est res] = lsqnonlin(@(guess)vsasl_fun...
    (guess, aqparms, parms, showFit, timeseries), ...
    guess0, ...
    lb,ub, opts);
%}
est

return

function res = vsasl_fun(guess, aqparms, parms, showFit, data )
% function res = VSASL_input_fun(del1, del2, del3, D, CBVa)


del1 = aqparms.del1;
del2 = aqparms.del2;
del3 = aqparms.del3;
ro_times = aqparms.ro_times;

f = guess(1);
CBVa = guess(2);
D = guess(3);


DURATION = 5;
t = linspace(0,DURATION,DURATION*1000);
dt = t(2)-t(1);

% Physiological constants
M0 = parms.M0;
alpha = parms.alpha;
R1a = 1/1.6;  % sec-1
R1t = 1/1.4;  % sec-1
lambda = 0.9;
M0 = parms.M0;
alpha = parms.alpha;
R1a = parms.R1a;
R1t = parms.R1t;  % sec-1
lambda = parms.lambda;


%Make an input function for VSASL
inp = zeros(size(t));
inp(round(del1/dt):round((del1+del2)/dt)) = 1;

decay = alpha*exp(-R1a*(t-del1));
decay(1:round(del1/dt)) = 0;
gkernel = t.* exp(-(t/D*2).^2);
gkernel = gkernel/sum(gkernel);
Sa = conv(inp, gkernel);
Sa = Sa(1:length(inp)).*decay;

% form a tissue retention and decay function
tis_fun = exp(-(R1t + f/lambda ) *(t));
tis_fun = tis_fun / max(tis_fun);

% The tissue signal is the convolution of the arterial input into the
% tissue. It assumes that the whole voxel feeding the tissue
St = f*conv( Sa , tis_fun);
St = St * dt; % adjust into units of seconds
St = St(1:length(inp));

% NOW We adjust the arterial signal to the arterial blood volume
Sa = Sa * CBVa;

signal = Sa+St;
ro_times = (del1+del2+del3+ro_times);
signal = signal(round(ro_times/dt));
res = signal-data;

if showFit
    hold off
    plot(t, Sa,'r','LineWidth',2);
    hold on
    plot(t, St, 'b','LineWidth',2);
    plot(t, St+Sa, 'g','LineWidth',2);
    plot(ro_times, signal, 'ko')

    if data~=0
        plot(ro_times, data, 'r*')
    end

    legend(...
        'Arterial Signal','Capillary/Tissue signal','Both', ...
        'Location','NorthEast')

    line([1 1]*(del1+del2), [0 max(Sa)])
    line([1 1]*(del1), [0 max(Sa)])
end

return

