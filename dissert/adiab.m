function result = adiab()

% Solves simultaneous system of diff. eq 
% Here the system of sequations will be the Bloch Equations with
% added terms to account for Magnetization transfer and for a
% variable frequencey RF sweep.

	% Define constants for the problem:   Milk !!	
	global Grad Vel T1 T2 H1 gamma M0 M0a M0b;
	global T1b  T2b ka kb ;

	Grad = 780; 	% mGauss/cm
	Vel = 20;	% cm/sec.
	gamma = 4.258;	% Hz/mGauss
	H1 = 100;	% mGauss (It takes 60 mGauss to satisfy
			% 	Gz*Vz <= gamma * H1 ^2
	ka =  1.02;	% Rate Ma -> Mb
	%f = 0.001;
	f =  0.02;	% Ratio Mb/Ma
	kb =  f*ka;	% Rate Mb -> Ma
	
	% First pool of protons:free
	T1 = 1.75;	% sec.
	T2 = 0.063;	% sec.
	%T2 = 0.163;	% sec.
	M0a = [0; 0; 1];	% initial condition: relaxed magnetization

	% Second pool of protons: bound
	T1b = 0.5;	% sec.
	%T2b = 0.001;	% sec.
	T2b = 0.00023;
	M0b = [0; 0; 1];	% initial condition: relaxed magnetization
		

	% Concatenating both vectors into one ...
	M0 = [	M0a ; 
		M0b];



	% Examine a range of velocities
	v = [1 5 10 15 20];
	for i=1:5
		Vel = v(i);
		iterate(strcat('vel', num2str(Vel)) );
	end;

	% eliminate MT effect
	% (removing these lines will bring it back)
	ka = 0;
	kb = 0;
	iterate('nomt');

	% use T1sat instead of MT effect in sec.  (T1sat of milk)
	T1 = 0.98;		
	iterate('t1sat');

return

%%%%%%%%%%%%%%%

function result =iterate(string1)
	
	output = diff_eq;
	write_mat(output, strcat(string1,'.dat') );

	hold on;
	grid on;
	plot(output(:,1),output(:,2),'b')
	plot(output(:,1),output(:,3),'r')
	title(string1);

return


%%%%%%%%%%%%%%%

function result = diff_eq()

% This is where the differential equation solver is called
% The file Bloch2.m contains the definition of the equations of two
% pools of protons with magnetization transfer
% The file Bloch.m contains the definition of the equations for
% a single pools of protons.

	global Grad Vel T1 T2 H1 gamma M0 M0a M0b;
	global T1b  T2b ka kb ;

%	options = odeset('RelTol',0.01);
%	[T,M2,S2] = ode23s('bloch2', [-1 1], M0, options);

	options = odeset('RelTol',0.01, 'InitialStep', 0.01);
	[T,M2,S2] = ode23s('bloch2', [-1:0.001:1], M0, options);
%	[T,M2,S2] = ode23s('bloch', [-1:0.01:1], M0a, options);
	result = [T(:,1)  M2(:,3)  M2(:,6)];
%	result = [T(:,1)  M2(:,3)];


return

%%%%%%%%%%%%%%%%%%%%%%

function result = not_used()

% This doesn't have ANY use.  
% It's just a place to cut and paste sections of code 
% that can be used under different circumstances

	Vel = 15;
	output = diff_eq('Vel = 15')
	write_mat(output, 'Vel15.dat');

	Vel = 10;
	output = diff_eq('Vel = 10')
	write_mat(output, 'Vel10.dat');

	Vel = 5;
	output = diff_eq('Vel = 5')
	write_mat(output, 'Vel05.dat');

	Vel = 1;
	output = diff_eq('Vel = 1')
	write_mat(output, 'Vel01.dat');



	output = diff_eq('using T1sat')
	write_mat(output, 't1sat.dat');
	
		figure;
	hold on;

	[T,M2,S2] = ode23s('bloch2', [-1:0.01:1], M0, options);
	result = [T(:,1)  M2(:,3)  M2(:,6)];
	write_mat(result, 'res2.dat');
	subplot(2,2,1),	plot(T,M2(:,3),'b')

	[T,M2,S2] = ode23s('bloch2', [-1:0.001:1], M0, options);
	subplot(2,2,2),plot(T,M2(:,3),'g')
	result = [T(:,1)  M2(:,3)  M2(:,6)];
	write_mat(result, 'res3.dat');

	[T,M2,S2] = ode23s('bloch2', [-1:0.0001:1], M0, options);
	subplot(2,2,3),plot(T,M2(:,3),'r')
	result = [T(:,1)  M2(:,3)  M2(:,6)];
	write_mat(result, 'res4.dat');

	[T,M2,S2] = ode23s('bloch2', [-1:0.00001:1], M0, options);
	subplot(2,2,4),plot(T,M2(:,3),'v')
	result = [T(:,1)  M2(:,3)  M2(:,6)];
	write_mat(result, 'res5.dat');


return
