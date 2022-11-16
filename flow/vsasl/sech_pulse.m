pl = 4; % ms
dt = 0.002; % ms
t = [-pl/2:dt:pl/2];
gambar = 4.258; %kHz/g
gam = 2*pi*gambar;

% hyperbolic secant pulse
mu = 5;
beta = 2;
a = 1;
amp = gam*a*sech(beta.*t);
freq = -mu*beta*tanh(beta.*t);

%{
beta=0.5;
amp = gam*a*cos(beta.*t) ;
freq = -mu*beta*sin(beta.*t); 
%}

ph = cumsum(freq).*dt;
ph = ph - ph(round(length(ph)/2));
rf = amp.*exp(i*ph);

BW = mu*beta/2/pi % in kHz

% plots
subplot(221)
plot(t,abs(rf))
subplot(222)
plot(t,freq/gam)
subplot(223)

plot(amp,freq)

beff = amp + i*freq/gam;
dthetadt = -diff(angle(beff))/dt;
dthetadt = [dthetadt 0];
kappa = abs(beff)./dthetadt;

subplot(224)
plot(t,kappa)
