deltaM = 0.01;
M0 = 1;
r1a = 1/1.67;
r1tis = 1/1.4;
label_dur = 3;
PLD = 1.6;
lambda = 0.9;
lambda = 1;
transit = 1.2;
TR = 5.0;
alpha = 0.9;

% PCASL solution from Buxton model according to white paper
f = deltaM  *lambda * r1tis * exp(PLD * r1a) ...
    / (2 * M0 * alpha * (1 - exp(-label_dur *r1a)))


% make the input function:
dt = 0.001;
inp = zeros(10/dt, 1);
m = zeros(size(inp), 1);
t = linspace(0, 10, 10/dt);

inp(transit/dt : (transit+label_dur)/dt) = 2*M0*alpha*exp(-transit*r1a);

% acquire images at these points:
aqs = [(label_dur+PLD) ,  TR+(label_dur+PLD)]./dt;
aqs = round(aqs);

% differential equation method (no MT included - should cancel out)
for n=2:length(inp)
    dm =  f * inp(n-1) - f*m(n-1)/lambda - m(n-1)*r1tis;       
    m(n) = m(n-1) + dm*dt;
end

% Sanity check:  Using convolution to calculate concentration of label in tissue:
ttmp = linspace(0,10, 10/dt)';
ret = exp( -(f/lambda + r1tis) * ttmp);
mt =  M0 * f * conv(inp , ret) * dt;

% clip the results
mt = mt(1:length(t));
inp = inp(1:length(t));

figure(13)
plot(t,mt); hold on
plot(t, m, 'g');
plot(t,inp,'r');
plot(t(aqs),mt(aqs),'g*');
hold off
legend('Conv.', 'Diff. Eq', 'Input Fun', 'AQ') 

new_deltaM  = -diff(mt(aqs))
