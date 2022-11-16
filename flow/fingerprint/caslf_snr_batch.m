function [BAT_error, flow_error, t1_error, cbv_error] = caslf_snr_batch(timing_parms)
Nreps = 20;
num_SNR_levels = 5;
SNRvals =  linspace(10,10^4,num_SNR_levels);

num_SNR_levels = 1;
SNRvals = 1000;


BAT = zeros(num_SNR_levels,Nreps);
flows = zeros(num_SNR_levels,Nreps);
cbvs = zeros(num_SNR_levels,Nreps);
T1s = zeros(num_SNR_levels,Nreps);

doFigs = 1;

%  generates dictionary

dict_phys_parms.r1tis = linspace(0.2,2,10);
dict_phys_parms.flip =   deg2rad(linspace(10,40,5));
dict_phys_parms.bat = linspace(0.5, 3, 10);
dict_phys_parms.bat2 = 0.3;
dict_phys_parms.f = linspace(0,100,20) / 6000;
dict_phys_parms.cbva = 0.01; % [linspace(0.005, 0.03, 5) ];
dict_phys_parms.kfor = 0.02; % [0 0.02 0.04 0.06];
dict_phys_parms.Disp = [40];
dict_phys_parms.mtis0 = 1;

fprintf('\nGenerating dictionary for SNR analysis ...')
[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

save flow_dictionary.mat dict parms

% generate an "observed"  signal 
my_phys_parms = parms(end/2-6);

% then alter them a bit...
 my_phys_parms.flip =   deg2rad(40);
 my_phys_parms.f = 0.005;
 my_phys_parms.bat = 1.2;
 my_phys_parms.r1tis = 1.3 ;
% my_phys_parms.bat2 = 0.3;
% my_phys_parms.cbva = 0.01;
% my_phys_parms.kfor = 0.02; 
% my_phys_parms.Disp = 40;
% my_phys_parms.mtis0 = 1;

% or select some existing set of paramaters from the main dictionary

dofigs = 0;
doSub = 0;
obssignal = gen_signals_150521(my_phys_parms , timing_parms, dofigs, doSub);

% 12/2/15 chop off 10 samples from begining - it screws everything up!
obssignal = obssignal(11:end);
dict = dict(:, 11:end);

for s=1:length(SNRvals)
    parfor n=1:Nreps
        
        fprintf('\nTesting different SNR levels ...SNR =  %f, rep n. = %d', SNRvals(s), n)
        
        % Adding noise to the observation
        NoiseLevel = abs(mean(obssignal))./SNRvals(s);
        tmpraw = obssignal + (NoiseLevel).*randn(size(obssignal));
        
        % searching for match
        best = flex_search_150521(dict, parms, tmpraw, 1,  [], timing_parms);
        
        BAT(s,n) = parms(best).bat;
        flows(s,n) = parms(best).f;
        cbvs(s,n) = parms(best).cbva;
        T1s(s,n) = 1/parms(best).r1tis;
        
    end
end

% figure out the errors in the match:
BAT_error = abs(100*(BAT - my_phys_parms.bat)/my_phys_parms.bat);
BAT_error = mean(BAT_error,2);

flow_error = abs(100*(flows - my_phys_parms.f)/my_phys_parms.f);
flow_error = mean(flow_error, 2);

cbv_error = abs(100*(cbvs - my_phys_parms.cbva)/my_phys_parms.cbva);
cbv_error = mean(cbv_error, 2);

t1_error = abs(100*(T1s - 1/my_phys_parms.r1tis)*my_phys_parms.r1tis);
t1_error = mean(t1_error, 2);

BAT_var = var(BAT,[],2);
flow_var = var(flows,[],2);
cbv_var = var(cbvs,[],2);
t1_var = var(T1s, [],2);

save SNR_errors.mat *error *var SNR*

%{
figure (26)

subplot(231), plot((SNRvals), (flow_error));
title('  MAE of flow')
xlabel('SNR')
ylabel('flow % MAE')
axis([0 max(SNRvals) -10 100])

subplot(232), plot((SNRvals), (BAT_error));
title('MAE of BAT Time (not log)')
xlabel('SNR')
ylabel('BAT Time % MAE')
axis([0 max(SNRvals) -10 100])

subplot(233), plot((SNRvals), (t1_error));
title(' MSE of T1 ')
xlabel('SNR')
ylabel('T1 % MSE')
axis([0 max(SNRvals) -10 100])


subplot(234), plot(SNRvals, (flow_var));
title('  Variance of Flow estimate')
xlabel('SNR')
ylabel('Flow variance')

subplot(235), plot(SNRvals, (BAT_var));
title('Variance of BAT Time estimate')
xlabel('SNR')
ylabel('BAT Time variance')

subplot(236), plot(SNRvals, (t1_var));
title(' Variance of CBV estimate')
xlabel('SNR')
ylabel('T1 variance')
%}
return

