% check that we're using the right scale factor for the jump size.
% calculate the random walk for a particle using the following 
% step size:

Nspins = 50;   % Number of particles to simulate
so = 100;
G=1000;   %Gauss/m
radius = 40e-9;    %(m)
T=0.05;     % time period to simulate (sec)
D= radius^2/(2*T)   % calculate the maximum allowed ADC.
ADC=D;
GAMMA = 4258*2*pi;	        %rad/s/G	

ADC = 0.94e-9    % (m^2/s)

dt = 1e-5;
Nsteps = T/dt;

dist = sqrt(2*ADC*dt*Nsteps);
scale = dist/sqrt(Nsteps);

%x=randn(Nsteps,Nspins) * scale/sqrt(2);
%y=randn(Nsteps,Nspins) * scale/sqrt(2);


% make a gradient
t = [1:Nsteps]*dt;
grad = (t >=  (Nsteps/4)*dt).*(t <= (Nsteps/2*dt)).* (-G) ... 
         + (t >= (Nsteps/2)*dt).*(t <= (3*Nsteps/4)*dt).* G ;

for c=1:Nspins
	fprintf('\rspin number %d',c);
	[x y] = bounce(Nsteps,radius,scale);
    phi = GAMMA * dt * cumsum(grad .* x);

end


endx = xx(end,:);
endy = yy(end,:);
dist2= sqrt(endx.^2 + endy.^2);


err = dist - dist2;

hist(err)
hist(dist2,50)
max(dist2)
M=[];
phi = cumsum(xx,1) ;
phi= sum(phi,2).*grad*dt;
%plot(phi)
magnet = ones(size(phi)).* exp((i* phi));
M = [M; magnet]; 

MM=sum(M);
subplot(211),plot(angle(MM))
subplot(212),plot(abs(MM))



% the maximum distance traveled by the particle is the radius
% of the circle.  The Apparent D is:
so = 100;
radius = 40e-9;    %(m)
T=0.1;
D= radius^2/(2*T)
s = so*exp(-bval*D)

% in vitro ADC for choline:
ADC = 0.94e-9    % (m^2/s)
s3=so*exp(-bval*ADC)

change=(s-s2)/so *100



% plot the change in signal with respect to ADC
so=100;
bval= 7e9; % s/m^2   =   7e3 s/mm^2

ADC = 1e-16 : 1e-14 : 1e-9; 
s = so*exp(-bval*ADC);
ds_dD = -bval*so*exp(-bval*ADC);

subplot(211),plot(ADC,s)
index=find( abs(ADC - 0.94e-9) < 1e-15)
hold on
plot(ADC(index), s(index),'or')
index=find( abs(ADC - D) < 1e-14)
plot(ADC(index), s(index),'or')
subplot(212),plot(ADC, ds_dD)
