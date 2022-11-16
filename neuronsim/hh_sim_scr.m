% How to use it
% =============
% 1) reset variables by
%    >> clear
% 2) set/override parameters
%    Vfix = 0 (normal mode) / 1 (keeps voltage at a fixed value)
%    V = capacitor voltage [mV] at time t=0
%    T = simulation time [ms]
%    dt = simulation step [ms]
%    t = initial time [ms]
%    Iext = external current [uA] at every time point (or at least in the interval
%           where it is non-zero)
% 3) run
%    >> hh
%
% Example
% >> clear
% >> hh
% >> Iext=pulse(0:dt:T,2,1,-10); t=0; hh
% >> Iext=0; % turn off the external current
% >> hh
% >> % reproduction of Fig.20 in Hodgkin,Huxley, J.Physiol., 500-544, (1952).
% >> V=-15; hh
%
if ~exist('cont')
	clear all
	% steady state conductances in m.mho/cm2
	gna0 = 120;
	gk0 = 36;
	gl = 0.3;
	Gconduct=0.0001; % 10MOhm
	
	% potentials are in milivolts
	Vna = -115;
	Vk = 12;
	Vl = -10.613;
	dVdt = 0;
	
	% membrane cap.
	Cm = 1.0; % 1.0 uF/cm2
		
	% initial values for gating constants:
	n = 0.86;  % potassium
	m = 0.9;  % sodium
	h = 0;    % sodium
	
	Iext=0; % no external current
	Vfix=0; % flowing potential
	dt=0.01; % msec - time step
	T=100;  % msec - duration
	V=0;
	t=0;
	
	cont=0;
    
    dt = 0.001;
    T = 50;
    
    
end

nsteps=ceil(T/dt);

% continue
t=t(end); 
t(nsteps+1)=0; % declare size
V=V(end);
V(nsteps+1)=0; % declare size
len=length(Iext);
if len<=nsteps
	Iext(nsteps+1)=0; % declare size
	Iext(len+1:nsteps+1)=0;
end

% use these for debugging:
h2 = zeros(1,nsteps);
m2 = h2;
n2 = m2;

for s=1:nsteps
	
	% V(s) = V(s) +  Vapplied(s);
	% Na channel conductance change rates:
	Am = 0.1 * (25+V(s)) / ( exp( (25+V(s))/10 )  - 1) ;
	Bm = 4 * exp(V(s)/18);

	Ah = 0.07 * exp(V(s)/20);
	Bh = 1/(exp( (30+V(s)) / 10) + 1);

	dmdt = Am * (1-m) - Bm * m;
	dhdt = Ah * (1-h) - Bh * h;

	h = h + dhdt * dt;
	m = m + dmdt * dt;

	gna = (m^3) * h * gna0;


	% K channel conductance change rates:
	An = 0.01*(10+V(s)) / ( exp( (10+V(s))/10 )  - 1) ;
	Bn = 0.125*exp(V(s)/80);

	dndt = An*(1-n) - Bn*n;
	n = n + dndt*dt;

	gk = (n^4) * gk0;
	
    % currents through the branches of the network
	Il  = gl  * (V(s)-Vl);
	Ik  = gk  * (V(s)-Vk );
	Ina = gna * (V(s)-Vna);
	Icond = Gconduct*(V(s));
    
    % current through capacitor:
    Ic= Iext(s) - ( Ina + Ik + Il + Icond);

    if Vfix
		V(s+1) = V(s);
	else		
		V(s+1) = V(s) + Ic*dt/Cm;
	end
	t(s+1) = t(s) + dt;
    
    % use these for debugging:
    h2(s) = h; m2(s) = m; n2(s) = n;
	
end

if 0
    disp(['V=',num2str(V(end)),'mV   dV/dt=',num2str(Ic*dt/Cm),'mV/ms   t=',num2str(t(s)),'ms'])
    disp(['Ic=',num2str(Ic),'uA   Ina=',num2str(Ina),'uA   Ik=',num2str(Ik),'uA'])
    disp(['                  Il=',num2str(Il),'uA   Icond=',num2str(Icond),'uA'])
    disp(['gna=',num2str(gna),'mS   gk=',num2str(gk),'mS'])
    disp(['Ah=',num2str(Ah),'1/ms   Bh=',num2str(Bh),'1/ms   h=',num2str(h),'   dh/dt=',num2str(dhdt),'1/ms'])
    disp(['Am=',num2str(Am),'1/ms   Bm=',num2str(Bm),'1/ms   m=',num2str(m),'   dm/dt=',num2str(dmdt),'1/ms'])
    disp(['An=',num2str(An),'1/ms   Bn=',num2str(Bn),'1/ms   n=',num2str(n),'   dn/dt=',num2str(dndt),'1/ms'])
end

if cont==0
	% set parameters for normal run
	T=50;
	dt=0.01;
    
    dt = 0.001;
    
	cont=1;
else
	subplot(212)
    plot(t,-V,'r'); title('mem. voltage [mV]')
    %axis([t(1) t(end) -100 150] )
    subplot(211)
    plot(t,-Iext,'b'); 
    %axis([t(1) t(end) -100 100] )
    xlabel('time [ms]'); title('-Iext [uA]');
end
