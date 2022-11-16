function result = Bloch(t, M)
% Defines the Bloch equations for use with ODE23 solution method
% there is no MTC anywhere

	% constants for the problem
	global Grad Vel T1 T2 H1 gamma M0;	

	result(1,1) = -M(1)/T2 + gamma*Grad*Vel*t*M(2);
	result(2,1) = -M(2)/T2 - gamma*Grad*Vel*t*M(1) + gamma * H1 * M(3);
	result(3,1) = -(M(3) - M0(3))/T1 - gamma* H1 * M(2);

return;
