% script to illustrate the properties of the different subtraction methods
% in ASL
figure

TR=2;
nyq=1/(2*TR);
tlen = 500;
NITER = 100;

simp = zeros(NITER, tlen/2);
run = zeros(NITER,tlen-1);
sinterp=zeros(NITER,tlen);
orig=zeros(NITER,tlen);


for iter=1:NITER
    fprintf('\riteration:  %d', iter')
    %noise = make1overf(tlen, nyq);
%    noise = rand(tlen,1);
    noise = makeAR1(tlen,0.6);
    raw = zeros(tlen,1);
    %raw=(1:tlen)';
    %raw(1:2:end)=2;
    raw = raw + noise ;
    %raw(1:2)=1;
    
    orig(iter,:) =  raw';
    t=0:tlen-1;
    
    % method 1: simple subtraction
    
    xt = raw(1:2:end);
    xc = raw(2:2:end);
    
    asl = xt - xc;
    simp(iter,:) =  asl';
    
    
    % method 2: running subtraction
    asl=zeros(tlen-1,1);
    for n=2:length(raw)
        asl(n-1) = (raw(n)-raw(n-1)) * (-1)^n;
    end
    run(iter,:)= asl';
    
    
    % method 3: sinterp subtraction
    % (A) using my own kernel
    xt = raw;
    xc = raw;
    
    xt(2:2:end) = 0;
    xc(1:2:end) = 0;
    
    xt=[xt(2:end); 0];
    
    kernel = 0.5* sinc((1.5/pi)*[-100:100]);
    
    xt2 = conv(xt,kernel);
    xc2 = conv(xc,kernel);
    xt2 = xt2(100:end-101);
    xc2 = xc2(100:end-101);
    asl = xt2 - xc2;
    
%     % (B) using Matlab's kernel:
%     xt = raw(1:2:end);
%     tt =t(1:2:end);
%     xc = raw(2:2:end);
%     tc = t(2:2:end);
%     
%     xt = interp1(tt,xt, t, 'sinc', 'extrap');
%     xc = interp1(tc,xc, t, 'sinc', 'extrap');
%     asl = xt' - xc';
%  
%     % (C) in Frequency domain:
%     xt = raw(1:2:end);
%     tt =t(1:2:end);
%     xc = raw(2:2:end);
%     tc = t(2:2:end);
%     w = [-pi:2*pi/length(xc):pi];
%     w=w(1:end-1)';
%     fxt = fftshift(fft(xt));
%     fxc = fftshift(fft(xc)).*exp(i*w);
%     fasl = fxt-fxc;
%     fasl = [zeros(length(fasl)/2,1); fasl ; zeros(length(fasl)/2,1)];
%     asl = ifft(fftshift(fasl));
%   
    sinterp(iter,:) = asl';
    
end

w= [-nyq:2*nyq/(length(raw)-1):nyq];
%w= [-2*nyq:2*nyq/(length(raw)-1):2*nyq];

% xo=[];
% for count=1:size(orig,1)
%     tmp=xcorr(orig(count,:));
%     xo = [xo ; tmp];
% end
% mforig = mean(fftshift(abs(fft(xo,[],2))),1);

mforig = mean(fftshift(abs(fft(orig,[],2))),1);

subplot(421), plot(w,mforig ), title('raw'),
axis tight




% simple subtraction:
nyq=1/(2*TR);
% 
% xs=[]
% for count=1:size(simp,1)
%     tmp=xcorr(simp(count,:));
%     xs = [xs ; tmp];
% end

%mfsimp = mean(fftshift(abs(fft(xs,[],2))),1);
%mfsimp = [zeros(1,tlen/2)  mfsimp   zeros(1,tlen/2)];
mfsimp = mean(fftshift(abs(fft(simp,[],2))),1);
mfsimp = [zeros(1,tlen/4)  mfsimp   zeros(1,tlen/4)];

f= [-nyq:2*nyq/(length(raw)-1):nyq];
%f= [-2*nyq:2*nyq/(length(raw)-1):2*nyq];

subplot(423)
plot (f, mfsimp), title('simple subtraction')
axis tight

subplot(424)
w = f*pi/nyq;
rect=zeros(size(w));
rect(length(rect)/4:3*length(rect)/4)=1;
plot(w, (mfsimp./mforig),'r');
hold on
% I don't know about that 0.5 factor....
plot(w, 0.5* abs((1-exp(-i*w/2)+ 1 - exp(-i*(w-2*pi)/2)) ).*rect)
axis ([-pi pi 0 5])




% running subtraction

nyq=1/(2*TR);

% xr=[]
% for count=1:size(run,1)
%     tmp=xcorr(run(count,:));
%     xr = [xr ; tmp];
% end
% mfrun = mean(fftshift(abs(fft(xr,[],2))),1);
% f= [-2*nyq:2*nyq/(length(raw)-1):2*nyq];
% subplot(425), plot (f(2:end-1),mfrun), title('running subtraction')
% subplot(426), plot(w(2:end-1), abs(mfrun./mforig(2:end-1)),'r');

mfrun = mean(fftshift(abs(fft(run,[],2))),1);
f= [-nyq:2*nyq/(length(mfrun)-1):nyq];
subplot(425), plot (f,mfrun), title('running subtraction')
axis tight

subplot(426), plot(w(1:end-1), abs(mfrun./mforig(1:end-1)),'r');
hold on
plot(w,abs(1-exp(-i.*(w-pi))))
axis([ -pi pi 0 5])






% sinc subtraction

% xi=[]
% for count=1:size(sinterp,1)
%     tmp=xcorr(sinterp(count,:));
%     xi = [xi ; tmp];
% end
% mfsinterp = mean(fftshift(abs(fft(xi,[],2))),1);

mfsinterp = mean(fftshift(abs(fft(sinterp,[],2))),1);
f= [-nyq : 2*nyq/(length(mfsinterp)-1) : nyq];

subplot(4,2,7)
plot (f, mfsinterp), title('sinc subtraction')
axis tight

subplot(4,2,8), 
w = f*pi/nyq;
plot(w, mfsinterp./mforig,'r');
hold on, 
subplot(4,2,8)
plot(w, 0.5* abs((1-exp(-i*w/2)- 1 + exp(-i*(w-2*pi)/2)) ).*rect)
axis([ -pi pi 0 5])




return

% testing analytical solutions 1/f noise
% simple subtraction:
w=-2*pi:0.1:2*pi

x1 = abs(2./w).*(1-exp(j*w/2));
x2 = abs(2./(w-2*pi)).*(1-exp(j*(w-2*pi)/2));

x2(end/2:end)=0;
x3 = abs(2./(w-2*pi)).*(1-exp(j*(w-2*pi)/2));
x3(1:end/2)=0;

plot(abs(x1))
plot(abs(x2))
plot(abs(x3))

x = x1+x2+x3;
rect=zeros(size(x));
rect(end/4:end*3/4)=1;
plot(w, abs(x).*rect)





