


%%%% now do the experimental data
figure;  set(gcf,'Position',[1 1 450,450]);
warning off

TR=1.4;
nyq=1/(2*TR);
%%%%%%
load ('voxels.mat');
tlen = size(mytimeseries_act,1);
NITER = size(mytimeseries_act,2);
types=2;

simp = zeros(NITER, tlen/2);
run = zeros(NITER,tlen-1);
sinterp=zeros(NITER,tlen);
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
X = [Xt Xc X];

Bhat_raw = [];
Bhat_simp = [];
Bhat_run = [];
%%%

for noisetype=1:types
    Bhat1=[];
    Bhat2=[];
    Bhat3=[];
    for iter=1:NITER
        fprintf('\rnoise type : %d , pixel:  %d', noisetype, iter')
        switch(noisetype)

            case 1
                title_str = 'Resting Voxels';
                % use data in the file: resting case
                noise = mytimeseries_rest(:,iter);
                t=my_glm(X,noise,[0 0 1]');
                load glm
                % remove the effect
                % noise= noise -beta_est(2:end)*X - beta_est(1);

            case 2
                title_str = 'Active Voxels';
                % use data in the file: active voxels case
                noise = mytimeseries_act(:,iter);
                t=my_glm(X,noise,[0 0 1]');
                load glm
                %  remove the effect from the data
                % noise= noise - beta_est(2:end)*X - beta_est(1) ;
        end

        t=0:tlen-1;

        % Process the data with all three methods and
        % do the GLM estimates

        % method 1: no preprocessing
        raw = noise - mean(noise);
        orig(iter,:) =  raw';

        D=D1;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat1 = [Bhat1; tmp];

        % method 2: simple subtraction

        xt = raw(1:2:end);
        xc = raw(2:2:end);

        asl = xt - xc;
        simp(iter,:) =  asl';

        D=D2;
        tmp = pinv(X'*D'*D*X)*X'*D'*D*raw;
        Bhat2 = [Bhat2; tmp];


        % method 3: running subtraction
        asl=zeros(tlen-1,1);
        for n=2:length(raw)
            asl(n-1) = (raw(n)-raw(n-1)) * (-1)^n;
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
        Bhat_sur = [Bhat_sur; tmp'];

    end

    % tabulate the Bhat for each data type:
    Bhat_raw = [Bhat_raw Bhat1];
    Bhat_simp = [Bhat_simp Bhat2];
    Bhat_run = [Bhat_run Bhat3];
    Bhat_sur = [Bhat_sur Bhat4];


    % Average the frequency responses of the processing
    % over pixels

    w = [-length(raw)/2:length(raw)/2-1]*2*pi/length(raw);
    f = w*nyq/(pi);

    mforig = mean(fftshift(abs(fft(orig,[],2))),1);
    %    nrg = sum(mforig);
    %    mforig = mforig/nrg;


    subplot(3,types, types*0+ noisetype)
    plot(f,mforig), title(title_str,'FontWeight','Bold','FontSize',18)
    axis([0 nyq 0 1.5*max(mforig(end/2+2:end)) ])



    % simple subtraction:

    mfsimp = mean(fftshift(abs(fft(simp,[],2))),1);
    %mfsimp = mfsimp/sum(mfsimp(end/2:end));
    %mfsimp = mfsimp/sum(mfsimp);
    mfsimp = [zeros(1,tlen/4)  mfsimp   zeros(1,ceil(tlen/4))];


    subplot(3,types,types*1 + noisetype)
    plot (f, mfsimp), title('Pairwise subtraction','FontWeight', 'Bold','FontSize',18)
    axis([0 nyq 0 1.5*max(mfsimp(end/1.8 +1:end)) ])



    % running subtraction

    w2 = [-length(run)/2:length(run)/2-1]*2*pi/length(run);
    f2 = w2*nyq/(pi);

    mfrun = mean(fftshift(abs(fft(run,[],2))),1);
    %    mfrun = mfrun/sum(mfrun(end/2:end));
    %    mfrun = mfrun/sum(mfrun);
    subplot(4,types ,types*2 + noisetype)
    plot (f2,mfrun), title('Running subtraction','FontWeight', 'Bold','FontSize',18)
    axis([0 nyq 0 1.5*max(mfrun(end/1.9:end)) ])
    fatlines


    % Surround subtraction

    mfsur = mean(fftshift(abs(fft(sur,[],2))),1);
    w2 = [ -length(mfsur)/2 : length(mfsur)/2-1 ]  * 2*pi/length(mfsur);
    f2 = w2*nyq/pi;
    subplot(4,types ,types*4 + noisetype)
    plot (f2,mfsur), title('Surround subtraction','FontWeight', 'Bold')
    fatlines


    % Now calculate the frequency responses base on the
    % analytical solutions

    % no Processing:
    x = mforig;

    subplot(4,types,noisetype)
    hold on
    %plot(f,abs(x),'r')
    legend('Original Data')
    legend boxoff
    fatlines

    %     % simple subtraction
    % y(w) = x(w/2) .* (1 - exp(-i*w/2))
    ysimp = x.*(1-exp(-i*w));
    % aliasing:  Note that it's very important that you apply the
    % 1-exp(-i*w) factor before you do the aloasing part...
    ysimp_alias = fftshift(ysimp);
    ysimp = (ysimp + ysimp_alias)/2;
    ysimp(1:end/4) = 0;
    ysimp(3*end/4:end) = 0;

    subplot(4,types,types+noisetype)
    hold on,
    plot(f, abs(ysimp),'r')
    fatlines
    legend('Numeric','Analytical')
    legend boxoff

    % running subtraction
    k = 1-exp(-i*(w+pi));
    yrun = k.*fftshift(x);

    subplot(4,types,types*2+noisetype)
    hold on
    plot(nyq*w/pi, abs(yrun),'r')
    fatlines
    legend('Numeric','Analytical')
    legend boxoff
end



%
%
%
% function noise = makeAR1(len, rho)
% % function noise = makeAR1(length, rho)
% %
% % returns a noise vector of 'length' smaples
% % with AR(1) characteristics, sampled at 'sampling_rate' (in Hz)
% % the vector is normalized such that the variance is one and the
% % mean is zero
%
%
%
% %
% % length = 1200;
% % rho=0.9
%
% noise = zeros(len,1);
% noise(1)=rand(1,1);
% for count=2:len
%     noise(count) = noise(count-1)*rho + rand(1,1);
% end
%
% noise = noise/count;
% % make sure that the variance is one and the mean is zero
% noise = noise - mean(noise);
% noise = noise / sqrt(var(noise));
%
% % subplot(211),plot(noise)
% % subplot(212), plot(abs(fft(xcorr(noise))))
%
% return