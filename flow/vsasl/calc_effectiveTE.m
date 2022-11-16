function effTE = calc_effectiveTE(basename, b1scale, Gmax)
%
%function effTE = calc_effectiveTE(basename, b1scale, Gmax)
% 
% units: Gauss
% simulate pulse on stationary tissue
% calculate how much signal is lost at the end of the pulse
% using an arbitrary T2
% calculate what TE would produce that amount of signal loss
%
str     = sprintf('%s.theta.txt', basename);
phs     = load(str);
str     = sprintf('%s.rho.txt', basename);
b1mag   = load(str);
str     = sprintf('%s.grad.txt', basename);
Gz      = load(str);

zpos = 0;

b1mag   = b1mag * b1scale/32766;
phs     = phs * pi/32766;
Gz     = Gz * Gmax/32766;

RF = b1mag .* exp(i*phs);

Mi = [0 0 1];

Bx = real(RF);
By = imag(RF);
Bz = Gz.*zpos;


nstep = length(RF);
dt = 4e-6; 

beff = [Bx By  Bz];

% unit conversionsl  Blochsim works in Tesla and mili seconds
T1 = 1.680;
T2 = 0.080;

T1 = T1 * 1e3;      % s to ms
T2 = T2 * 1e3;      % s to ms
dt = dt * 1e3;      % s to ms
beff = beff*1e-4; % G to Tesla

M = blochsim(Mi, beff, T1, T2, dt, nstep);

M = M(1:length(RF),:);

Mpost = M(end,:)

effTE = -T2*log(abs(Mpost(3)))


t = [0:length(Gz)-1] * dt;

subplot(411)
area(t, abs(RF)*1000);
grid off
ylabel('B1 amplitude (mG)')
axis tight
set(gca, 'box','off')

subplot(412)
plot(t, angle(RF),'k'); 
grid off
ylabel('B1 phase (rad)')
axis tight
set(gca, 'box','off')

subplot(413)
area(t, Gz)
grid off
ylabel('Grad.(G/cm)')
axis tight
set(gca, 'box','off')

subplot(414)
plot(t, M(:,3),'k'); 
grid off
ylabel('M_z')
axis tight
set(gca, 'box','off')
