	% Define constants for the problem:   Milk !!	
	global Grad Vel T1 T2 H1 gamma M0 M0a M0b;
	global T1b  T2b ka kb ;

	Grad = 780; 	% mGauss/cm
	Vel = 20;	% cm/sec.
	gamma = 4.258;	% Hz/mGauss
	H1 = 100;	% mGauss (It takes 60 mGauss to satisfy
			% 	Gz*Vz <= gamma * H1 ^2
	T1 = 1.5;	% sec.
	T2 = 0.053;	% sec.
	M0 = [0; 0; 1];	% initial condition: relaxed magnetization

	v = [1 5 10 15 20 30];
	allM = [];
	for i=1:6
		Vel = v(i)
		options = odeset('RelTol',0.01);
		[T,M2,S2] = ode23s('bloch', [-1 :0.0001: 2], M0(1:3), options);
		plot(T, M2(:,3)); drawnow ; hold on
		allM = [allM M2(:,3)];
		whos M2
	end;
	figure
	plot(T, allM); drawnow
	fatlines
	dofontsize(16)
	xlabel('Time (sec)')
	ylabel('Z Magnetization')
	legend('1' ,'5' ,'10', '15', '20', '30')
	axis([-1 2 -1.1 1.1])
