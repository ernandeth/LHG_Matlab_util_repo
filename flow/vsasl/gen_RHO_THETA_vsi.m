dacmax = hex2dec('7ffe');

NPOINTS = 50;
rf1 = hanning(NPOINTS);
tmp = abs(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
plot(tmp)


% first write it as intt16
filename = sprintf('myHanning%d.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);


% then use GE software to put in proper format
str = sprintf('!xlatebin -o myHanning%d.rho myHanning%d.bin', NPOINTS, NPOINTS);
eval(str)


%%
NPOINTS = 50;

rf1 = ones(NPOINTS,1);
tmp = abs(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);
plot(tmp)


% first write it as intt16
filename = sprintf('myRect%d.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);


% then use GE software to put in proper format
str = sprintf('!xlatebin -o myRect%d.rho myRect%d.bin', NPOINTS, NPOINTS);
eval(str)


%%
%Now a hyperbolic secant
dacmax = hex2dec('7ffe');

Nsegs = 1;
segN = 50;
NPOINTS  = segN *Nsegs;
dt = 4e-6;
t = linspace(0, 50e-3, NPOINTS); % ms;
Amp = 50e-7 ; % Tesla  (50 mGauss)

gambar = 42576e3;     % gamma in Hz/T
sweepWidth = 500; % Hz

beta = 1;
k = 0.8;

weights = sech(beta * linspace(-pi , pi, NPOINTS));
weights = weights - min(weights);
weights = weights / (max(weights) - min(weights));

% frequency sweep of the pulse train envelope
RF_dw = tanh(k*beta * linspace(-pi , pi, NPOINTS)); % works well
RF_dw = RF_dw * 2/( max(RF_dw)-min(RF_dw));

B_eff = complex(gambar*weights*Amp,   sweepWidth*RF_dw);

prec_rate = gambar*abs(B_eff) ;   % in Hz
prec_rate = prec_rate(1:end-1);
tip_rate = diff(angle(B_eff)) / dt;
adiabty = prec_rate ./ tip_rate;


rf1 = weights .* exp(i* pi*RF_dw*sweepWidth .* (t - t(end/2)) );


min(adiabty)

subplot(231)
plot(weights)
title('|B1| weights')

subplot(232)
plot(RF_dw)
title('dw (Hz)')

subplot(233)
plot(angle(rf1))
title('RF phase (rads)')

subplot(234)
plot(adiabty)
title('Adiabaticity')

subplot(235)
plot(B_eff)
title('B_{eff} trajectory')
axis square

out= rf1;

subplot(236)
plot(abs(out))
hold on
plot(angle(out), 'r')
title('RF pulse')
hold off

NPOINTS = length(out);
% Now write out scanner files:
tmp = abs(out);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);

% first write it as intt16
filename = sprintf('mysech_%d.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o mysech_%d.rho mysech_%d.bin', NPOINTS, NPOINTS);
eval(str)

% Now write out scanner files:
tmp = angle(out);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);

% first write it as intt16
filename = sprintf('mysech_%d.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o mysech_%d.theta mysechS_%d.bin', NPOINTS, NPOINTS);
eval(str)


%%  the VSI pulses
dacmax = hex2dec('7ffe');

Nsegs = 1;
segN = 50;
NPOINTS  = segN *Nsegs;
dt = 4e-6;
t = linspace(0, 50e-3, NPOINTS); % ms;
Amp = 50e-7 ; % Tesla  (50 mGauss)

gambar = 42576e3;     % gamma in Hz/T
sweepWidth = 500; % Hz

beta = 1;
k = 0.8;

weights = sech(beta * linspace(-pi , pi, NPOINTS));
weights = weights - min(weights);
weights = weights / (max(weights) - min(weights));

% frequency sweep of the pulse train envelope
RF_dw = tanh(k*beta * linspace(-pi , pi, NPOINTS)); % works well
RF_dw = RF_dw * 2/( max(RF_dw)-min(RF_dw));

B_eff = complex(gambar*weights*Amp,   sweepWidth*RF_dw);

prec_rate = gambar*abs(B_eff) ;   % in Hz
prec_rate = prec_rate(1:end-1);
tip_rate = diff(angle(B_eff)) / dt;
adiabty = prec_rate ./ tip_rate;


rf1 = weights .* exp(i* pi*RF_dw*sweepWidth .* (t - t(end/2)) );


min(adiabty)

subplot(231)
plot(weights)
title('|B1| weights')

subplot(232)
plot(RF_dw)
title('dw (Hz)')

subplot(233)
plot(angle(rf1))
title('RF phase (rads)')

subplot(234)
plot(adiabty)
title('Adiabaticity')

subplot(235)
plot(B_eff)
title('B_{eff} trajectory')
axis square

out = [];
zseg = zeros(1,segN);
for n=1:Nsegs
    beg=(n-1)*segN +1;
    fin = segN*n;
    out = [out rf1(beg:fin) zseg ];
end
rf1=out(1:end-segN);

% The gradient pulses
out = [];
Gseg = zeros(1,segN);
legN = floor((segN-2)/4);
Gseg(2 : legN+1) = linspace(0,1,legN);
Gseg(legN+2 : 2*legN+1) = linspace(1,0,legN);
Gseg( 2*legN+2 : 3*legN+1) = linspace(0,-1,legN);
Gseg(3*legN+2 : 4*legN+1) = linspace(-1,0,legN);

G = zeros(size(rf1));
for n=1:Nsegs
    beg=(n-1)*segN +1;
    fin = segN*n;
    out = [out G(beg:fin) Gseg ];
end
G=out(1:end-segN)';


subplot(236)
area(G)
hold on
plot(abs(rf1))
plot(angle(rf1), 'r')
title('RF pulse')
hold off

NPOINTS = length(rf1);

% Now write out scanner files:
tmp = abs(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);

% first write it as intt16
filename = sprintf('myVSI_%d.rho.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.rho myVSI_%d.rho.bin', NPOINTS, NPOINTS);
eval(str)

% Now write out scanner files:
tmp = angle(rf1);
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);


% first write it as intt16
filename = sprintf('myVSI_%d.theta.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.theta myVSI_%d.theta.bin', NPOINTS, NPOINTS);
eval(str)


% write out a text file with the gradient waveform
str = sprintf('save myVSIGrad_%d.txt G -ascii', NPOINTS)
eval(str)

% Now write out scanner files:
tmp = G ;
tmp = tmp * dacmax / max(tmp);
tmp = 2*round(tmp/2);

% first write it as intt16
filename = sprintf('myVSI_%d.grad.bin',NPOINTS);
fp = fopen(filename,'wb','b');
fwrite(fp, tmp, 'int16');
fclose(fp);

% then use GE software to put in proper format
str = sprintf('!xlatebin -o myVSI_%d.grad myVSI_%d.grad.bin', NPOINTS, NPOINTS);
eval(str)

