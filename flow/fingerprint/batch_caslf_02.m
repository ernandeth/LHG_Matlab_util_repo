
Nreps = 50;
SNRvals =  linspace(10,500,10);
trans = zeros(10,Nreps);
flows = zeros(10,Nreps);
cbvs = zeros(10,Nreps);

% first run: generates dictionary

[dictionary parms] = gen_dictionary

myparms = parms(1);
myparms.f = 0.012;
myparms.transit = 1.2;

obssignal = gen_signals(myparms, 1 , 1);

useFileDictionary = 1;

for s=1:length(SNRvals)
    for n=1:Nreps
       
       % Adding noise to the observation
        NoiseLevel = abs(mean(obssignal))./SNR;
        obssignal = obssignals + (NoiseLevel).*randn(size(obssignals));

        bestparms = search_dictionary(obssignal, useFileDictionary)
        
        trans(s,n) = bestparms.transit;
        flows(s,n) = bestpartms.f;
		cbvs(s,n) = bestparms.cvba;
    end
end

trans_error = (100*(trans - 1.2)/1.2).^2;
trans_error = mean(trans_error,2);
flow_error = (100*(flows - 0.01)/0.01).^2;
flow_error = mean(flow_error, 2);
cbv_error = (100*(cbvs - 0.05)/0.05).^2;
cbv_error = mean(cbv_error, 2);

trans_var = var(trans,[],2);
flow_var = var(flows,[],2);
cbv_var = var(cbvs,[],2);

save milf_sim_0524

figure(6)
subplot(231), plot(SNRvals, (flow_error));
title('  MSE of flow')
xlabel('SNR')
ylabel('flow MSE')

subplot(232), plot(SNRvals, (trans_error));
title('MSE of Trans Time (not log)')
xlabel('SNR')
ylabel('Trans Time MSE')

subplot(233), plot(SNRvals, (cbv_error));
title(' MSE of CBV ')
xlabel('SNR')
ylabel('CBV MSE')


subplot(234), plot(SNRvals, (flow_var));
title('  Variance of Flow estimate')
xlabel('SNR')
ylabel('Flow variance')

subplot(235), plot(SNRvals, (trans_var));
title('Variance of Trans Time estimate')
xlabel('SNR')
ylabel('Trans Time variance')

subplot(236), plot(SNRvals, (cbv_var));
title(' Variance of CBV estimate')
xlabel('SNR')
ylabel('CBV variance')


