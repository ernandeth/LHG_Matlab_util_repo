function result= quick_mrf_sensitivity(timing_parms)
% function result= quick_mrf_sensitivity(timing_parms)
%
% Now test the schedule's sensitivity by varying one parameter at a tim
% we cut them in half.
%

parms.mtis0 =     1 ;
parms.Disp =      40;
parms.r1blood = 1/1.7;

parms.f =       0.01;
parms.cbva =    0.01;
parms.bat =     0.1;
parms.r1tis =   1/0.9;%1.4;
parms.flip =    deg2rad(40);
parms.r2tis=    1/0.090;

doSub = 0;

for n=1:5
    parms2 = parms;
    switch(n)
        case 1
            str = 'CBF'
            parms2.f = parms.f/2;
        case 2
            str = 'CBV'
            parms2.cbva = parms.cbva/2;
        case 3
            str = 'BAT'
            parms2.bat = parms.bat/2;
        case 4
            str = 'R_1'
            parms2.r1tis = parms.r1tis/2;
        case 5
            str = 'R_2'
            parms2.r2tis = parms.r2tis/2;
    end
    
    test_signal = abs(single(  ...
        gen_signals_vs_220421(parms,...
        timing_parms,...
        0, doSub))); % ...
    
    test_signal2 = abs(single(  ...
        gen_signals_vs_220421(parms2,...
        timing_parms,...
        0,doSub))); ...
        figure(2)
    
    
    rms = mean(abs((test_signal-test_signal2)./(test_signal))) ;
    
    doSub = 0;
  
    subplot(5,2,2*n-1)
    plot(test_signal); hold on;
    plot(test_signal2); hold off
    title(str)
    
    subplot(5,2,2*n)
    plot( (test_signal-test_signal2)/norm(test_signal) )
    title(sprintf('NRMS change  : %0.2e',...
        rms));
    
    switch(n)
        case 1
            result.df = rms;        
        case 2
            result.dcbva = rms;
        case 3
            result.bat = rms;        
        case 4
            result.r1tis = rms;
        case 5
            result.r2tis = rms;
    end
    
end 
    