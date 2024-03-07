function result= quick_mrf_sensitivity(timing_parms, parms)
% function result= quick_mrf_sensitivity(timing_parms, parms)
%
% Now test the schedule's sensitivity by varying one parameter at a tim
% we cut them in half.
%


doSub = 0;

for n=1:5
    parms2 = parms;
    switch(n)
        case 1
            str = 'CBF'
            parms2.f = parms.f*2;
        case 2
            str = 'CBV'
            parms2.cbva = parms.cbva*2;
        case 3
            str = 'BAT'
            parms2.bat = parms.bat*2;
        case 4
            str = 'R_1'
            parms2.r1tis = parms.r1tis/2;
        case 5
            str = 'R_2'
            parms2.r2tis = parms.r2tis/2;
    end
    
    
        %gen_signals_vs_220421(parms, ...
        %gen_signals_vs_230321(parms, ...

    test_signal = abs(single(...
        gen_signals_vs_230918(parms, ...
        timing_parms, ...
        0, doSub, 1e-3)));

        %gen_signals_vs_220421(parms2,...
        %gen_signals_vs_230321(parms2,...
    test_signal2 = abs(single(...
        gen_signals_vs_230918(parms2, ...
        timing_parms,...
        0,doSub, 1e-3))); 
    
        figure(2)
    
    
    nrms = norm(test_signal-test_signal2) / mean(test_signal) ;
    
  
    subplot(5,3,3*n-2)
    plot(test_signal); hold on;
    plot(test_signal2); hold off
    title(str)
    
    subplot(5,3,3*n-1)
    plot( (test_signal-test_signal2) )
    title(sprintf('NRMS change  : %0.2e',...
        nrms));

    subplot(5,3,3*n)
    plot(test_signal(1:2:end)-test_signal(2:2:end))
    hold on
    plot(test_signal2(1:2:end)-test_signal2(2:2:end) )
    title('Pairwise diffs (time)')
    hold off

   subplot(5,3,3*n)
   plot(abs(fftshift(fft(test_signal-test_signal2))));
   title('FFT of change')

    switch(n)
        case 1
            result.dSdf = nrms;        
        case 2
            result.dSdcbva = nrms;
        case 3
            result.dSdbat = nrms;        
        case 4
            result.dSdr1tis = nrms;
        case 5
            result.dSdr2tis = nrms;
    end
   
    total_duration = ...
        sum(timing_parms.del1(:)) ...
        + sum(timing_parms.del2(:)) ...
        + sum(timing_parms.del3(:))...
        + sum(timing_parms.RO_time(:))
end 
   
[total_duration result.dSdf result.dSdcbva result.dSdbat result.dSdr1tis result.dSdr2tis] 

