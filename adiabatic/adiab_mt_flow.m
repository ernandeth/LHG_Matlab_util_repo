% Solves simultaneous system of diff. eq 
% Here the system of sequations will be the Bloch Equations with
% added terms to account for Magnetization transfer and for a
% variable frequencey RF sweep.
global parms

	Grad = 780; 	% mGauss/cm
	Vel = 11.2;	% cm/sec. Experimental
	gamma = 4.258;	% Hz/mGauss
	H1 = 100;	% mGauss (It takes 60 mGauss to satisfy
			% 	Gz*Vz <= gamma * H1 ^2
	water = 1
	
    T1sat = 1.36;	
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
	T2b = 280e-6;
        
    % Values from Stanisz' paper:
    f = 0.028;
    ka = 35;
    kb = ka*f;
    T2b = 280e-6;

% Water values: 
if water
	T1 = 2.3
	T2 = 0.365
	ka = 0
	kb=0
	Vel = 10
end
%%%%%%%%%%%%%%%%%
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
    
    R1 = 1/T1;
 
    parms =	[ Grad; Vel; T1; T2; H1; gamma; T1b;  T2b; ka; kb; M0];

    % Examine a range of velocities
    
        string1 = sprintf('MTVel%02f.dat', Vel);
        options = [];%odeset('RelTol',0.01);

        [T,M2,S2] = ode45('bloch_mtc', [-0.5 1], M0);
        output = [T(:,1)  M2(:,3)];

        hold on; grid on;
        plot(output(:,1),output(:,2),'k')
        title(string1);
        drawnow
        str = sprintf('save %s output -ASCII', string1);
        eval(str)
        
        alpha = (1 - min(M2(:,3)))/2;

    title('Measured T1 Decay During Transit');
    xlabel('Time (s)'), ylabel('Mz/Mz0')
    dofontsize(16)
    fatlines
    %print -dtiff figure2B

   data = load('BloodPhantomData.dat')
   if water
   	data = load('WaterPhantomData.dat')
   end

   t = data(:,2);
   sig = data(:,1);
   plot(t,sig,'ok')
   hold on

   sigT1 = 1 - (1-sig(1))*exp(-(t-t(1))/T1);
   sigT1sat = 1 - (1-sig(1))*exp(-(t-t(1))/T1sat);
   if ~water
	plot(t,sigT1sat,'.-k')
   end
   plot(t,sigT1,'-*k')
	
   legend('Full Simulation','Measured Data','T1sat decay','T1 decay')
   if water
	legend('Full Simulation','Measured Data','T1 decay')
   end
   fatlines
   dofontsize(16)
   axis([0 0.5 -0.6 0])

mseT1 = sum((sig - sigT1).^2)
mseT1sat = sum((sig - sigT1sat).^2)
simt = output(:,1);
for n=1:length(t)
	suspects = find(abs(simt-t(n)) < 0.0001 );
	inds(n) = suspects(ceil(end/2));
end
tt = output(inds,1);
simSig = output(inds,2);
mseSim = sum((sig - simSig).^2)

return 
    % eliminate MT effect
    % (removing these lines will bring it back)
    ka = 0;
    kb = 0;
    parms =	[ Grad; Vel; T1; T2; H1; gamma; T1b;  T2b; ka; kb; M0];


    string1 = 'nomt';

    [T,M2,S2] = ode45('bloch_mtc', [-1 1], M0, options);
    output = [T(:,1)  M2(:,3)];

    hold on;     grid on;
    plot(output(:,1),output(:,2),'--k')
    drawnow
    str = sprintf('save %s output -ASCII', string1);
    eval(str)


