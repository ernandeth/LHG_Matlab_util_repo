function result = bloch_mtc(t, M)

% Defines the Bloch equations for use with ODE23 solution
% (Uses a Runge-Kutta approximation method)
% This version includes two pools of exchanging protons.
global parms

	% constants for the problem
    M0 = zeros(6,1);

    Grad = parms(1);
    Vel = parms(2);
    T1 = parms(3);
    T2 = parms(4);
    H1 = parms(5);
    gamma = parms(6);
	T1b  = parms(7);
    T2b = parms(8);
    ka = parms(9);
    kb = parms(10);
    M0(:)= parms(11:16);
    
	% Free pool:
    result = zeros(6,1);
    
	result(1) = -M(1)/T2 + gamma*Grad*Vel*t*M(2)  ;			% Mx(t)
	result(2) = -M(2)/T2 - gamma*Grad*Vel*t*M(1) + gamma * H1 * M(3);	% My
	result(3) = -(M(3) - M0(3))/T1 - gamma*H1*M(2) + ka*M(6) - kb*M(3);	% Mz

	% Bound pool :

	result(4) = -M(4)/T2b + gamma*Grad*Vel*t*M(5)   ;			% Mx(t)
	result(5) = -M(5)/T2b - gamma*Grad*Vel*t*M(4) + gamma * H1 * M(6);	% My
	result(6) = -(M(6) - M0(6))/T1b - gamma*H1*M(5) - ka*M(6) + kb*M(3);	% Mz
    


return;
