function result= quick_mrf_model_test(timing_parms)
% function result= quick_mrf_model_test(timing_parms)
%
% compares two gen_signals models ata time
%

parms.Mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01;
parms.cbva =    0.01;
parms.bat =     0.1;
parms.r1tis =   1/0.9;%1.4;
parms.flip =    deg2rad(40);
parms.r2tis=    1/0.090;
parms.b1err = 0;

doSub = 0;
%{
for n=1:5
    parms2 = parms;
    switch(n)
        case 1
            str = 'CBF'
            parms2.f = parms.f/22;
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
    
%}  
        %gen_signals_vs_220421(parms, ...
        %gen_signals_vs_230321(parms, ...

    test_signal = abs(single(...
        gen_signals_vs_230718(parms, ...
        timing_parms, ...
        0, doSub, 1e-3)));

        %gen_signals_vs_220421(parms2,...
        %gen_signals_vs_230321(parms2,...
    test_signal2 = abs(single(...
        gen_signals_vs_230918(parms, ...
        timing_parms,...
        0,doSub, 1e-3))); 
    
        figure(2)
    
    
    rms = mean(abs((test_signal-test_signal2)./(test_signal))) ;
    
      
    subplot(2,1,1)
    plot(test_signal); hold on;
    plot(test_signal2); hold off
    title('raw signals')

    subplot(2,1,2)
    plot( (test_signal-test_signal2)/norm(test_signal) )
    title(sprintf('NRMS change  : %0.2e',...
        rms));
   
    total_duration = ...
        sum(timing_parms.del1(:)) ...
        + sum(timing_parms.del2(:)) ...
        + sum(timing_parms.del3(:))...
        + sum(timing_parms.RO_time(:))

    result = [test_signal; test_signal2];
end 
    