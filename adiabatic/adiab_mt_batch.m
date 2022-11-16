% Solves simultaneous system of diff. eq 
% Here the system of sequations will be the Bloch Equations with
% added terms to account for Magnetization transfer and for a
% variable frequencey RF sweep.
global parms

	Grad = 780; 	% mGauss/cm
	Vel = 10;	% cm/sec.
	gamma = 4.258;	% Hz/mGauss
	H1 = 100;	% mGauss (It takes 60 mGauss to satisfy
			% 	Gz*Vz <= gamma * H1 ^2
	
    T1sat = 1.01;	
	ka =  1/T1sat;	% Rate Ma -> Mb
	f =  0.1;	% Ratio Mb/Ma, also ka/kb
	kb =  ka*f;	% Rate Mb -> Ma
    
	
	% First pool of protons:free
	T1 = 1.75;	% sec.
	T1 = 2.0;	% sec.
	T2 = 0.098;	% sec.
   
    % Use measurements from the UM 7T lab:
    T1 = 1.84;
    T2 = 0.062;
    
	% Second pool of protons: bound
	T1b = 0.5;	% sec.
	T2b = 0.001;
        
    % Values from Stanisz' paper:
    f = 0.028;
    ka = 35;
    kb = ka*f;
    T2b = 280e-6;
%     T1 = 1.932;
%     T2 = 0.275;
%     T1b = 0.05;
%     T2b = 0.280;
%     
   % Concatenating both vectors into one ...
    M0a = [0; 0; 1-f];	% initial condition: relaxed magnetization
    M0b = [0; 0; f];	% initial condition: relaxed magnetization
    M0 = [	M0a ;
        M0b];
    
    % aside:  this is what the MT spectrum should look like:
    R1 = 1/T1;
    %experimental data from Blood at 7 T (NKI)
    dw = [15000
        12000
        10000
        8000
        5000
        4000
        2000
        1000
        500
        100];

    mzD = [0.977255105
        0.977772034
        0.960713363
        0.965882657
        0.927371414
        0.932282243
        0.828896356
        0.632721633
        0.457741018
        0.159989661];

    % Fit My version of Henkelmann's 1993 equation:
    % to get he ka and f paramaters
    MTparms = [H1 T1 T2 f T1b T2b];
    guess0 = [ka];
    [estimate , resnorm, res, ex, output]  = lsqnonlin(@MTspec_lsq,...
        guess0, [], [], [], dw,MTparms, mzD);
    
    mz = MTspec_lsq(estimate, dw,MTparms);
    
    plot(dw,mz, 'k')
    hold on, plot(dw,mzD,'ko')
    
    % from the literature?
    mz = MTspec_lsq([35], dw,MTparms);
    hold on, plot(dw,mz,'k')
    
    title('Magnetization Transfer Spectrum')
    legend(sprintf('fitted k = %2.2f',estimate), 'Blood Data');%, 'Stanisz k=35') 
    xlabel('Frequency (hz)'), ylabel('Mz/Mz0')
    dofontsize(16)
    fatlines

    % Here comes the simulation of the flowing spins:
    figure
    
    % Use this to consider the case without MT
    % ka = 0; kb=0;
    
    parms =	[ Grad; Vel; T1; T2; H1; gamma; T1b;  T2b; ka; kb; M0];

    % Examine a range of velocities
    v = [1 5 10 15 20 30];
    %v=10;
    %v=11.2;
    alphas = zeros(size(v));
    
    for i=1:length(v)
        parms(2) = v(i);
        Vel = v(i);
        string1 = sprintf('MTVel%02d.dat', Vel);
        options = [];%odeset('RelTol',0.01);

        [T,M2,S2] = ode45('bloch_mtc', [-1 1], M0);
        output = [T(:,1)  M2(:,3)];

        hold on; grid on;
        subplot 122, plot(output(:,1),output(:,2),'k')
        title(string1);
        drawnow
        str = sprintf('save %s output -ASCII', string1);
        eval(str)
        
        alphas(i) = (1 - min(M2(:,3)))/2;
    end;
    title('B. Effect of Velocity (MT present)');
    xlabel('Time (s)'), ylabel('Mz/Mz0')
    dofontsize(16)
    fatlines
    %print -dtiff figure2B

    % eliminate MT effect
    % (removing these lines will bring it back)
    ka = 0;
    kb = 0;
    Vel = 15 ; 
	parms =	[ Grad; Vel; T1; T2; H1; gamma; T1b;  T2b; ka; kb; M0];

    string1 = 'nomt';

    [T,M2,S2] = ode45('bloch_mtc', [-1 1], M0, options);
    output = [T(:,1)  M2(:,3)];

    hold on;     grid on;
    %plot(output(:,1),output(:,2),'--k')
    drawnow
    str = sprintf('save %s output -ASCII', string1);
    eval(str)


    % Make some figures from the simulations
    %figure
    mt = load ('MTVel20.dat');
    grid on;
    subplot 121, plot(mt(:,1),mt(:,2),'k')
    title('A. Effect of Magnetization Transfer');
    
    nomt = load ('NOMTVel20.dat');
    hold on;     grid on;
    plot(nomt(:,1),nomt(:,2),'--k')
    legend('MT present' ,'No MT')
    xlabel('Time(s.)')
    ylabel('Mz (a.u)')
    fatlines
    dofontsize(16)
    
    a1 = min(nomt(:,2));
    a1 = (1 - a1)/2;
    
    a2 = min(mt(:,2));
    a2 = (1 - a2)/2;
    
    a1end = (1-nomt(end,2))/2;
    a2end = (1-mt(end,2))/2;

    reduction = (a1 - a2)/a1 
    reduction2 = ((a1-a1end) - (a2-a2end)) / (a1-a1end)
    
    % make the plots based on the files
    for i=1:length(v)
        Vel = v(i);
        string1 = sprintf('MTVel%02d.dat', Vel);
        output = load(string1);

        hold on; grid on;
        subplot 122, plot(output(:,1),output(:,2),'k')
        title(string1);
        drawnow
        str = sprintf('save %s output -ASCII', string1);
        eval(str)

        alphas(i) = (1 - min(M2(:,3)))/2;
    end;
    title('B. Effect of Velocity (MT present)');
    xlabel('Time (s)'), ylabel('Mz/Mz0')
    dofontsize(16)
    fatlines
   print -dtiff figure01
   
   return 
    
    % Fit T1sat from progressive saturation experiment at NKI
    TR = [0.2 0.5 1 2 5 10];
    data = [1315.66 3266.67 5843.73 8825.4 11309.6 11469.1];

    Mo_guess = 1.2*data(end);
    T1_guess = 1;
    FlipFactor_guess = 1;

    LB = [100, 0.5 0 ];
    UB = [10e8, 3  2];
    optvar = optimset('lsqnonlin');

    guess0 = [Mo_guess; T1_guess; FlipFactor_guess];
    guess = lsqnonlin('sr_lsq', ...
        guess0, LB, UB, ...
        optvar, ...
        TR, ...
        data);

    M0= guess(1);
    T1sat = guess(2);
    FlipFactor = guess(3);


    % use T1sat instead of MT effect in sec.  
    T1 = T1sat;
    string1 = 't1sat';

    options = odeset('RelTol',0.01);

    [T,M2,S2] = ode45('bloch_mtc', [-1 1], M0, options);
    output = [T(:,1)  M2(:,3)+M2(:,6)];

    hold on;     grid on;
    plot(output(:,1),output(:,2),'k')
    title(string1);
    drawnow
    str = sprintf('save %s output -ASCII', string1);
    eval(str)
