close all
w = [-pi:0.01:pi-0.01];
nyq = 1;
tlen = length(w);
rho=0.99
% x = ones(size(w)) * mean(mforig);
x = ones(size(w));

% defining AR1 as:  x[n] = rho*x[n-1] + eps[n]
x = x ./ (1 - rho * exp(-j.*w));
x = x.^2;
% its derivative:
% FT{diff(y[n])} = w.* Y(w);
wx = w.*x ./ (1 - rho * exp(-j.*w));

% defining AR1 as:  eps[n] = rho*eps[n-1] +eps[n]
%x = x .* (1 + rho * exp(-j.*w));


%  w = [-length(x)/2:length(x)/2-1]*2*pi/length(x);
f = w*nyq/pi;


% identity

subplot(5,1,1)
hold on
plot(nyq*w/pi,abs(x),'r')
%plot(nyq*w/pi,abs(wx),'b')
plot(nyq*w/pi, abs((1-rho^2)./(1-2*rho*cos(w)+rho^2)),'g');
legend('data','derivative')
legend boxoff
title('AR1 and its derivative')
axis([0 nyq 0 2*abs(mean(x))])

% pairwise subtraction

ysimp = x.*(1-exp(-i*w));
% aliasing: Note that it's very important that you apply the
% 1-exp(-i*w) factor before you do the aliasing part...
ysimp = (ysimp + fftshift(ysimp) )/2;
ysimp(1:end/4) = 0;
ysimp(3*end/4:end) = 0;
ysimp = ysimp.^2;

%    ysimp = x + x(end:-1:1);
subplot(5,1,2)
hold on,
%plot(nyq*w/pi, abs(ysimp),'r')
ww=w(end/4+1:3*end/4-1);
yy=ysimp(end/4+1:3*end/4-1);
plot(nyq*ww/pi, abs(yy),'r')
title('pairwise')

legend('data','theoretical')
legend boxoff
axis([0 nyq 0 2*abs(mean(x))])



% Running subtraction
k = (1-exp(-i*(w+pi)))/2;
yrun = k.*fftshift(x);
yrun = yrun.^2;
subplot(5,1,3)
hold on
subplot(5,1,3)
plot(nyq*w/pi, abs(yrun),'r')
title('running')
legend('data','theoretical')
legend boxoff
axis([0 nyq 0 2*abs(mean(x))])


% surround subtraction:
k = (2 - exp(-i*(w+pi)) - exp(i*(w+pi)))/4;
ysur = k.*fftshift(x);
ysur = ysur.^2;
subplot(5,1,4)
hold on
plot(nyq*w/pi, abs(ysur),'r')
title('surround')
legend('data','theoretical')
legend boxoff
axis([0 nyq 0 2*abs(mean(x))])

% sinc subtraction
ysinc = x.*(1-exp(-i*w));
% aliasing: Note that it's very important that you apply the
% 1-exp(-i*w) factor before you do the aliasing part...
ysinc = (ysinc + fftshift(ysinc) )/2;
%ysinc(1:end/4) = 0;
%ysinc(3*end/4:end) = 0;

klen = 41;
t = linspace(-round(tlen/2),round(tlen/2),tlen);
s = sinc(t/2);
s(1:end/2-20)=0;
s(end/2+20:end) = 0;
skernel = fftshift(fft(s));
ysinc = ysinc .* skernel;
ysinc = ysinc.^2;
subplot(5,1,5)
hold on
plot(nyq*w/pi, abs(ysinc),'r')
title('sinc')

