close all; clear variables;

Nreps = 1000;
SNRvals =  [0.1 1:10:100];
trans = zeros(length(SNRvals),Nreps);
flows = zeros(length(SNRvals),Nreps);
T1s = zeros(length(SNRvals),Nreps);

% cbvs = zeros(length(SNRvals),Nreps);

doFigs = 0;

% first run: generates dictionary
load optimalTiming_100im_SNR75_2f.mat;
[dictionary parms] = gen_dictionary_140328(timing_parms);

myparms = parms(2710);
obssignal = gen_signals_140627(myparms, timing_parms , 0,0);

dict=abs(dict);
save dictionary dict parms ;


useFileDictionary = 1;

matlabpool('open',6)
for s=1:length(SNRvals)
    parfor n=1:Nreps
        % Adding noise to the observation
        NoiseLevel = mean(abs(obssignal))./SNRvals(s);
        fprintf('\rSNR = %f  Noise Level = %f',SNRvals(s),  NoiseLevel);
        
        tmp = abs(obssignal + (NoiseLevel).*randn(size(obssignal)));
        
        [bestparms] = search_dictionary(tmp, dict, parms, doFigs);
        
        trans(s,n) = bestparms.transit;
        flows(s,n) = bestparms.f;
        T1s(s,n) = 1./bestparms.r1tis;
    end
end
matlabpool('close')

trans_error = (100*abs(trans - myparms.transit)/myparms.transit);
trans_error2 = mean(trans_error,2);
flow_error = (100*abs(flows - myparms.f)/myparms.f);
flow_error2 = mean(flow_error, 2);
T1_error = (100*abs(T1s - 1./myparms.r1tis)/(1./myparms.r1tis));
T1_error2 = mean(T1_error, 2);

trans_var = std(trans_error,[],2);
flow_var = std(flow_error,[],2);
T1_var = std(T1_error,[],2);

% save milf_sim1

figure;
subplot(131), errorbar(SNRvals, (flow_error2),flow_var);
title('Flow: Mean Percent Error')
xlabel('SNR')
ylabel('Flow: Mean Percent Error')

subplot(132), errorbar(SNRvals, (trans_error2),trans_var);
title('Trans Time: Mean Percent Error')
xlabel('SNR')
ylabel('Trans Time: Mean Percent Error')

subplot(133), errorbar(SNRvals, (T1_error2),T1_var);
title('T1: Mean Percent Error')
xlabel('SNR')
ylabel('T1: Mean Percent Error')

save Timing_100im_SNR75_2f_sim trans flows T1s;
% temp=[flow_error, flow_var, trans_error, trans_var];
% SNRa = 5;
% NoiseLevel = mean(abs(obssignal))./SNRa;
% tmp = abs(obssignal + (NoiseLevel).*randn(size(obssignal)));
% [bestparms,entry] = search_dictionary(tmp', useFileDictionary, doFigs);
% figure;plot(1:25,abs(obssignal),1:25,tmp,1:25,dict(entry,:));legend('Original Signal','Signal,SNR=5','Matched Ditionary Entry');ylabel('Signal');xlabel('Images')

