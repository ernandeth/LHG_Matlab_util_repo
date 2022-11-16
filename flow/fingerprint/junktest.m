close all

y = load('gm_timecourse.dat');
y = differencer(y,3);
y = y-mean(y);
y = y/norm(y);

for beta=linspace(20,90,10);
    phys_parms.f =         0.01;
    phys_parms.mtis0 =     1 ;
    phys_parms.cbva =      0.03 ;
    phys_parms.transit=    0.2 ;
    phys_parms.kfor =     0.002 ;
    phys_parms.r1tis =     1/1.2  ;
    phys_parms.beta =      beta*pi/180 ; % flip angle in radians
    phys_parms.L = 1;
    phys_parms.Disp =      20;
    phys_parms.Ptime =     0.5;
    
    
    % Just a quick test to make sure everything works as expected
    % obs = gen_signals_150116(phys_parms, timing_parms, 1, doSub);
    
    obs = gen_signals_140618(phys_parms, timing_parms, 0, 0);
   obs = differencer(obs,3);
    obs = obs-mean(obs);
    obs = obs/norm(obs);
    
    hold on
    plot(obs)
    pause(0.1)
    drawnow
end

plot(y,'r');