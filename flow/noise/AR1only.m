% script to illustrate the properties of the different subtraction methods
% in ASL
close all;
warning off

rho=0.96;
TR=1.4;
nyq=1/(2*TR);
%%%%%%
load ('voxels.mat');
tlen =258;
NITER = 50;
types=1;

simp = zeros(NITER, tlen/2);
run = zeros(NITER,tlen-1);
sur = zeros(NITER,tlen-2);
sincc=zeros(NITER,tlen);
orig=zeros(NITER,tlen);

% differencing matrices
D1 = eye(tlen);
D2 = zeros(tlen/2,tlen);
for count=1:tlen/2
    D2(count,count*2-1)=1;
    D2(count,count*2)=-1;
end
D3 = zeros(tlen);
for count=1:tlen-1
    D3(count,count)=(-1)^(count-1);
    D3(count,count+1)=(-1)^(count);
end
D4 = zeros(tlen);
for count=1:tlen-2
    D4(count,count)=(-1)^(count-1);
    D4(count,count+1)=2*(-1)^(count);
    D4(count,count+2)=(-1)^(count-1);
end

% make up the sinc-interpolation matrix.
klen = 41;
t = linspace(-round(klen/2),round(klen/2),klen);
skernel = (sinc(t-0.5));
sk = zeros(1,length(skernel)*2);
sk(2:2:end) = skernel;
D5 = zeros(tlen + 2*klen);
for row=2:2:tlen+klen
    D5(row,row-1:row-2 + 2*klen ) = sk;

end
for row=3:2:tlen+klen
    D5(row,row-1:row-2 + 2*klen ) = -sk;
end

D5 = D5(klen:end-klen-1, 2*klen:end-klen+1);
for row=1:tlen
    % normalize the kernel
    D5(row,:) = D5(row,:) / abs(sum(D5(row,:)));
    % stick in the ones to represent the non-interpolated samples
    D5(row,row) = -(-1)^row;
end
D5 = D5(1:tlen, 1:tlen);

% design matrix
X = zeros(tlen,1);
X(round([18:20:360]/1.4))=1;
X = conv(X,spm_hrf(TR));
X = X(1:tlen);
% modulation of the effect
mm = ones(size(X));
mm(1:2:end)=-1;
X=X.*mm;
% baseline perfusion: regressor for controls and regressor for tag:
Xc = zeros(tlen,1);
Xt = zeros(tlen,1);
Xc(1:2:end)=1;
Xt(2:2:end)=1;
% put it all together:
X = [Xc Xt X];
%figure;  set(gcf,'Position',[10 10 650,450]);

imagesc(X)
colormap(gray);

Bhat_raw = [];
Bhat_simp = [];
Bhat_run = [];
Bhat_sur = [];
Bhat_sinc = [];

%%%
for noisetype=1:types
    %for noisetype=2:2
    Bhat1=[];
    Bhat2=[];
    Bhat3=[];
    Bhat4=[];
    Bhat5=[];

    for iter=1:NITER
        %fprintf('\rnoise type : %d , pixel:  %d', noisetype, iter')
        noise = makeAR1(tlen,rho);
        %                 noise = randn(tlen,1);
        %                 V = spm_Q(rho,tlen) ;
        % %                 % This is just a safety measure,to force homogeneous variance:
        %                  V = diag(1./sqrt(diag(V)))*V;
        %                  noise = V^1/2 *noise;


        raw = noise;% - mean(noise);
        orig(iter,:) =  raw';
        t=0:tlen-1;

        % method 1: no preprocessing
        D=D1;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat1 = [Bhat1; tmp];

        % method 2: pairwise subtraction

        xt = raw(1:2:end);
        xc = raw(2:2:end);

        asl = xt - xc;
        simp(iter,:) =  asl';

        D=D2;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat2 = [Bhat2; tmp];


        % method 3: Running subtraction
        asl=zeros(tlen-1,1);
        for n=2:length(raw)
            asl(n-1) = (raw(n)-raw(n-1)) / 2* (-1)^n;
        end
        run(iter,:)= asl';

        D=D3;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat3 = [Bhat3; tmp];

        % method 4: surround subtraction
        asl=zeros(tlen-2,1);
        for n=2:length(raw)-1
            asl(n-1) = (-raw(n+1)+ 2*raw(n) -raw(n-1)) / 4 * (-1)^n;
        end
        sur(iter,:)= asl';

        D=D4;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat4 = [Bhat4; tmp];

        % method 5: sinc subtraction
        %asl = sincsub(raw);
        asl = D5*raw;
        sincc(iter,:)= asl';

        D=D5;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat5 = [Bhat5; tmp];


    end

    % tabulate the Bhat for each data type:
    Bhat_raw = [Bhat_raw Bhat1];
    Bhat_simp = [Bhat_simp Bhat2];
    Bhat_run = [Bhat_run Bhat3];
    Bhat_sur = [Bhat_sur Bhat4];
    Bhat_sinc = [Bhat_sinc Bhat5];

    w = [-length(raw)/2:length(raw)/2-1]*2*pi/length(raw);
    f = w*nyq/pi;

    mforig = mean(fftshift(abs(fft(orig,[],2))),1);
%    mforig = fftshift(fft(mean(orig,1),[],2));
    
    subplot(5,types, types*0+ noisetype)
    plot(f,abs(mforig)), title(title_str,'FontWeight','Bold')
    axis([0 nyq 0 2*max(abs(mforig(end/2+15:end))) ])



    % pairwise subtraction:

    mfsimp = mean(fftshift(abs(fft(simp,[],2))),1);
%    mfsimp = mean(fftshift((fft(simp,[],2))),1);
    mfsimp = [zeros(1,ceil(tlen/4))  mfsimp   zeros(1,(tlen/4))];

    subplot(5,types,types*1 + noisetype)
    plot (f, mfsimp), title('Pairwise sub.','FontWeight', 'Bold')
    axis([0 nyq 0 2*max(abs(mfsimp(end/2 +15:end)) )])



    % Running subtraction
    w = [-length(run)/2:length(run)/2-1]*2*pi/length(run);
    f = w*nyq/pi;

    mfrun = mean(fftshift(abs(fft(run,[],2))),1);

    subplot(5,types ,types*2 + noisetype)
    plot (f,mfrun), title('Running sub.','FontWeight', 'Bold')
    axis([0 nyq 0 2*max(mfrun(end/2+15:end)) ])
    fatlines

    % Surround subtraction
    w = [-length(sur)/2:length(sur)/2-1]*2*pi/length(sur);
    f = w*nyq/pi;

    mfsur = mean(fftshift(abs(fft(sur,[],2))),1);

    subplot(5,types ,types*3 + noisetype)
    plot (f,mfsur), title('Surround sub.','FontWeight', 'Bold')
    axis([0 nyq 0 2*max(mfrun(end/2+15:end)) ])
    fatlines

    % Sinc subtraction
    w = [-length(sincc)/2:length(sincc)/2-1]*2*pi/length(sincc);
    f = w*nyq/pi;

    mfsinc = mean(fftshift(abs(fft(sincc,[],2))),1);

    subplot(5, types ,types*4 + noisetype)
    plot (f,mfsinc), title('Sinc sub.','FontWeight', 'Bold')
    axis([0 nyq 0 2*max(mfsinc(end/2+15:end)) ])
    fatlines


    %     Now do the analytical solutions
    %

    % theoretical normalized noise:
   
            x = ones(size(w));
            x = x ./ (1 - rho * exp(-j.*w));
            %x(ceil(length(x)/2)+1) = 0; % remove DC
            %x = x*abs(sum(mforig))/abs(sum(x)); % normalize to the energy in the simulation
            x = abs(x*sum(abs(mforig))/sum(abs(x))); % normalize to the energy in the simulation
            %x(130)=0;
            %x = mforig;%(end/2:end);


    %  w = [-length(x)/2:length(x)/2-1]*2*pi/length(x);
    f = w*nyq/pi;


    % identity

    subplot(5,types,noisetype)
    hold on
    plot(nyq*w/pi,abs(x),'r')
    %legend('data','theoretical')
    legend boxoff
    fatlines
    %axis([0 nyq 0 2*abs(mean(x))])

    % pairwise subtraction

    ysimp = x.*(1-exp(-i*w));
    % aliasing: Note that it's very important that you apply the
    % 1-exp(-i*w) factor before you do the aliasing part...
    ysimp = (ysimp + fftshift(ysimp) )/2;
    ysimp(1:end/4) = 0;
    ysimp(3*end/4:end) = 0;


    subplot(5,types,types+noisetype)
    hold on,
    plot(nyq*w/pi, abs(ysimp),'r')
    fatlines
    %legend('data','theoretical')
    %legend boxoff
    %axis([0 nyq 0 2*abs(mean(x))])

    plot(nyq*w/pi,mfsimp)


    % Running subtraction
    k = (1-exp(-i*(w+pi)))/2;
    yrun = k.*fftshift(x);
    subplot(5,types,types*2+noisetype)
    hold on
    plot(nyq*w/pi, abs(yrun),'r')
    fatlines
    %legend('data','theoretical')
    %legend boxoff
    %axis([0 nyq 0 2*abs(mean(x))])



    % surround subtraction:
    k = (2 - exp(-i*(w+pi)) - exp(i*(w+pi)))/4;
    ysur = k.*fftshift(x);
    subplot(5,types,types*3+noisetype)
    hold on
    plot(nyq*w/pi, abs(ysur),'r')
    fatlines
    %legend('data','theoretical')
    %legend boxoff
    %axis([0 nyq 0 2*abs(mean(x))])

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
    subplot(5,types,types*4+noisetype)
    hold on
    plot(nyq*w/pi, abs(ysinc),'r')
    fatlines
end


% tabulate the efficiencies into a table
mean_Bhat = [mean(Bhat_raw,1)  ;
    mean(Bhat_simp,1) ;
    mean(Bhat_run,1) ;
    mean(Bhat_sur,1) ;
    ];


mean_Bhat_rel = [mean(Bhat_raw,1) ./ mean(Bhat_simp,1)  ;
    mean(Bhat_simp,1)./ mean(Bhat_simp,1) ;
    mean(Bhat_run,1) ./ mean(Bhat_simp,1);
    mean(Bhat_sur,1) ./ mean(Bhat_simp,1);
    ];

std_Bhat = [std(Bhat_raw,1) ;
    std(Bhat_simp,1) ;
    std(Bhat_run,1) ;
    std(Bhat_sur,1) ;
    ];

std_Bhat_rel = [std(Bhat_raw,1) ./ std(Bhat_simp,1) ;
    std(Bhat_simp,1) ./ std(Bhat_simp,1) ;
    std(Bhat_run,1) ./ std(Bhat_simp,1) ;
    std(Bhat_sur,1) ./ std(Bhat_simp,1) ;
    ];