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

for noisetype=1:3
    
    for iter=1:NITER
        fprintf('\rNoise type : %d , iteration:  %d', noisetype, iter')
        switch(noisetype)
            case 1
                noise = rand(tlen,1);
                noise = noise -mean(noise);
            case 2
                noise = make1overf(tlen, nyq);
            case 3
                noise = makeAR1(tlen,0.9);
        end
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
    
    mforig = mean(fftshift(abs(fft(orig,[],2))),1);
    subplot(4,3, noisetype)
    plot(mforig), title('noise input')
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
    
    subplot(4,3,3 + noisetype)
    plot (f, mfsimp), title('simple subtraction')
    axis tight
    
    
    
    
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
    subplot(4,3 ,6 + noisetype)
    plot (f,mfrun), title('running subtraction')
    axis tight
    
    
    
    
    
    
    % sinc subtraction
    
    % xi=[]
    % for count=1:size(sinterp,1)
    %     tmp=xcorr(sinterp(count,:));
    %     xi = [xi ; tmp];
    % end
    % mfsinterp = mean(fftshift(abs(fft(xi,[],2))),1);
    
    mfsinterp = mean(fftshift(abs(fft(sinterp,[],2))),1);
    f= [-nyq : 2*nyq/(length(mfsinterp)-1) : nyq];
    subplot(4,3, 9 + noisetype)
    plot (f, mfsinterp), title('sinc subtraction')
    axis tight
    
   
end

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





