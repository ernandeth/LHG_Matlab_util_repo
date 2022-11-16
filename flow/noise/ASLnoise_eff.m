% script to illustrate the properties of the different subtraction methods
% in ASL - turboCASL in event related experiments.
%
% uses simulated data with AR(1) noise
%
% 4 differencing matrices.
% calculates power, efficiency, MSE, bias of the beta estimate as a
% function of SNR
%
batchmode=1;

if ~batchmode
    clear

    % set up the options here:
    DesignType = 2;  % 1=blocked, 2=ER Fixed, 3=ER randomized
    ResultsFile = 'ERFixed_sim_results';
    doSpectra = 0;
    showPlots = 0;
    doMeanCenter = 1;
    warning off
    doGLS = 0;  % other wise it's OLS
    rho = 0.3;
    TR = 1.4;
    nyq = 1/(2*TR);
    tlen = 258;
    NITER = 1;
    doNewMatrix = 1;  % generate a new design matrix or read it from a file?

    % determine a range of noise levels - SNRs..
    nlevels = [0.5 1:15 20 30];
    %nlevels = [0.5:0.5: 4];
    %nlevels = 5;
end

% Make the differencing matrices
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

%
Spower_SNR = zeros(5,length(nlevels));
eff_SNR = zeros(5,length(nlevels));
bias_SNR = zeros(5,length(nlevels));
std_SNR = zeros(5,length(nlevels));
bcrit_SNR = zeros(5,length(nlevels));
close all;

for nl=1:length(nlevels)
    fprintf('\rNoise level ...  %f ', nlevels(nl));

    noise_amp = nlevels(nl);

    simp = zeros(NITER, tlen/2);
    run = zeros(NITER,tlen-1);
    sinterp=zeros(NITER,tlen);
    orig=zeros(NITER,tlen);
    sincc=zeros(NITER,tlen);

    if doMeanCenter
        numBetas=3;
    else
        numBetas=2;
    end
    
    Bhat_raw = zeros(NITER,numBetas);
    Bhat_simp = zeros(NITER,numBetas);
    Bhat_run = zeros(NITER,numBetas);
    Bhat_sur = zeros(NITER,numBetas);
    Bhat_sinc = zeros(NITER,numBetas);

    VB_raw = zeros(NITER,numBetas);
    VB_simp = zeros(NITER,numBetas);
    VB_run = zeros(NITER,numBetas);
    VB_sur = zeros(NITER,numBetas);
    VB_sinc = zeros(NITER,numBetas);

    Spower_MC = zeros(5,NITER);
    eff_MC = zeros(5,NITER);

    %%%
    %%% do the estimation for each case

    for iter=1:NITER
        fprintf('\rSimulation number: %d noise = %f ', iter, nlevels(nl));


        if (doNewMatrix==1)
            fprintf(' Making new design matrix ...')
            % Make the design matrix (main effect)
            X = zeros(tlen,1);
            %X(round([18:10:360]/1.4))=1;

            
            % 1=blocked, 2=ER Fixed, 3=ER randomized
            switch DesignType
                case 1
                    isi = [1:30 60:90 120:150 180:210 240:tlen];
                case 2
                    isi = [1:20:360 2:20:360]/1.4;
                case 3
                    isi = rand(30,1)*12 + 5;
                    isi = round(cumsum(isi));
                    isi = isi(find(isi<tlen));
            end
            X(round(isi))=1;


            X = conv(X,spm_hrf(TR));
            X = X(1:tlen);
            X = X/max(X);

            % modulation of the effect
            mm = ones(size(X));
            mm(1:2:end)=-1;
            X=X.*mm/2;

            % baseline perfusion: regressor for controls and regressor for tag:
            Xc = zeros(tlen,1);
            Xt = zeros(tlen,1);
            Xc(2:2:end)=1;
            Xt(1:2:end)=1;

            % put in BOLD effect regressor here if desired

            % baseline signal
            Xb = ones(size(Xc));

            % put it all together:
            % note that in the baseline case we use only the tag to prevent
            % too much colinearity with the activation vector
            %X = [Xb Xc X];
            if doMeanCenter
                X = [Xc X Xb];
                X = X - repmat(mean(X), tlen, 1);
                B = [-1 -.5 100]';
            else
                X = [Xc X];
                B = [-1 -.5 ]';
            end

        else
            %    load /Users/hernan/matlab/flow/noise/ASLX.mat
            fprintf(' Reading design matrix ...')
            load /Users/hernan/matlab/flow/noise/ASLX2.mat
        end

        signal = X*B ;
        
        %title_str = 'AR(1) simulation';

        %noise = makeAR1(tlen,rho);
        %if doNewMatrix==1
        noise = randn(tlen,1);
        %end

        % a better way to create AR1 noise?
        % note - now we vary SNR by changing the amplitude of the effect
        V = spm_Q(rho,tlen) * noise_amp;
        % This is just a safety measure,to force homogeneous variance:
        V = diag(1./sqrt(diag(V)))*V;
        % W = V^(-1/2) to make inversions easier down the line...
        W = WKfun('mkW',[],V);

        %         D5W = zeros(size(W));
        %         for col=1:size(W,2)
        %             D5W(:,col) = sincsub(W(:,col));
        %         end
        %         D5V = zeros(size(V));
        %         for col=1:size(V,2)
        %             D5V(:,col) = sincsub(V(:,col));
        %         end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mean center the signal and noise cov V matrix here
        if doMeanCenter
            I = eye(size(V));
            J = ones(size(V))/length(signal);
            noise = V^(1/2) *noise;
            raw = noise + signal;
            raw    = (I-J)*raw;
            V = (I-J)*V*(I-J)';
        else
            noise = V^(1/2) *noise;
            raw = noise + signal;
        end

        orig(iter,:) =  raw';
        t=0:tlen-1;

        %         hold on, plot(D3*raw)
        %         drawnow

        %         D5raw = sincsub(raw);

        % method 1: no differencing
        D=D1;

        DX = D*X; 
        %DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_raw (iter,:)= tmp';

        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_raw = pinv(D*W*X)*D*W* V *W * D'*(pinv(D*W*X))';
            VB_raw(iter,:) = (diag(gls_VB_raw))';
        else
            ols_VB_raw = pinv(D*X)*D* V *D'*(pinv(D*X))';
            VB_raw(iter,:) = (diag(ols_VB_raw))';
        end

        % method 2: simple subtraction

        xt = raw(1:2:end);
        xc = raw(2:2:end);

        asl = xt - xc;
        simp(iter,:) =  asl';

        D=D2;
        
        DX = D*X; 
        %DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_simp (iter,:)= tmp';

        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_simp = pinv(D*W*X)*D*W* V *W * D'*(pinv(D*W*X))';
            VB_simp(iter,:) = (diag(gls_VB_simp))';
        else
            ols_VB_simp = pinv(D*X)*D* V *D'*(pinv(D*X))';
            VB_simp(iter,:) = (diag(ols_VB_simp))';
        end

        % method 3: running subtraction
        asl=zeros(tlen-1,1);
        for n=2:length(raw)
            asl(n-1) = (raw(n)-raw(n-1)) / 2 *(-1)^n;
        end
        run(iter,:)= asl';

        D=D3;
        
        DX = D*X; 
        %DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_run(iter,:)= tmp';

        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_run = pinv(D*W*X)*D*W* V *W * D'*(pinv(D*W*X))';
            VB_run(iter,:) = (diag(gls_VB_run))';
        else
            ols_VB_run = pinv(D*X)*D* V *D'*(pinv(D*X))';
            VB_run(iter,:) = (diag(ols_VB_run))';
        end

        % method 4: surround subtraction
        asl=zeros(tlen-2,1);
        for n=2:length(raw)-1
            asl(n-1) = (-raw(n+1)+ 2*raw(n) -raw(n-1)) / 4 * (-1)^n;
        end
        sur(iter,:)= asl';

        D=D4;

        DX = D*X; 
        %DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_sur(iter,:)= tmp';

        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_sur = pinv(D*W*X)*D*W* V *W * D'*(pinv(D*W*X))';
            VB_sur(iter,:) = (diag(gls_VB_sur))';
        else
            ols_VB_sur = pinv(D*X)*D* V *D'*(pinv(D*X))';
            VB_sur(iter,:) = (diag(ols_VB_sur))';
        end

        % method 5: sinc subtraction
        D=D5;
        
        DX = D*X; 
        %DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;

        Bhat_sinc(iter,:)= tmp';

        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_sinc = pinv(D*W*X)*D*W* V *W * D'*(pinv(D*W*X))';
            VB_sinc(iter,:) = (diag(gls_VB_sinc))';
        else
            ols_VB_sinc = pinv(D*X)*D* V *D'*(pinv(D*X))';
            VB_sinc(iter,:) = (diag(ols_VB_sinc))';
        end

        %         tmp = pinv(D5X'*D5X)*D5X'*D5raw;
        %         Bhat_sinc(iter,:)= tmp';
        %
        %         ols_VB_sinc = pinv(D5X)* D5V *(pinv(D5X))';
        %         % recall that W = V^(-1/2)
        %         gls_VB_sinc = pinv(D5W*X)*D5W * V *D5W'*(pinv(D5W*X))';





        % tabulate the stats for Monte Carlo
        mean_Bhat = ...
            [mean(Bhat_raw,1)  ;
            mean(Bhat_simp,1) ;
            mean(Bhat_run,1) ;
            mean(Bhat_sur,1) ;
            mean(Bhat_sinc,1) ;
            ]


        mean_Bhat_rel = mean_Bhat/mean(Bhat_simp,1);

        %     Monte Carlo
        %   std_Bhat = ...
        %         [std(Bhat_raw,0,1) ;
        %         std(Bhat_simp,0,1) ;
        %         std(Bhat_run,0,1) ;
        %         std(Bhat_sur,0,1) ;
        %         std(Bhat_sinc,0,1) ;
        %         ];
        %

        if doGLS==1
            % now using the GLS formula instead of the brute force ...
            std_Bhat = ...
                [diag(gls_VB_raw)';
                diag(gls_VB_simp)';
                diag(gls_VB_run)';
                diag(gls_VB_sur)';
                diag(gls_VB_sinc)';]
        else
            std_Bhat = ...
                [diag(ols_VB_raw)';
                diag(ols_VB_simp)';
                diag(ols_VB_run)';
                diag(ols_VB_sur)';
                diag(ols_VB_sinc)';
                ]

        end
        std_Bhat_rel = std_Bhat /std(Bhat_simp,0,1);

        eff = 1./std_Bhat;
        eff_rel = eff / eff(2);

        % bias and MSE on the third regressor (activation regressor)
        %     bias = [ Bhat_raw(:,3) Bhat_simp(:,3) Bhat_run(:,3) Bhat_sur(:,3) Bhat_sinc(:,3) ] ;
        %     bias = (bias - B(3)) / B(3);
        %     mean_bias = mean((bias),1)
        %     mean_bias_rel = mean_bias / mean_bias(2)
        %
        %     MSE = [ Bhat_raw(:,3)  Bhat_simp(:,3)  Bhat_run(:,3) Bhat_sur(:,3) Bhat_sinc(:,3)] ;
        %     MSE = (MSE - B(3)).^2 / (B(3))^2;
        %     mean_MSE = mean((MSE),1);

        bias = [ Bhat_raw(:,2) Bhat_simp(:,2) Bhat_run(:,2) Bhat_sur(:,2) Bhat_sinc(:,2) ] ;
        bias = (bias - B(2)) / B(2);
        mean_bias = mean((bias),1);
        mean_bias_rel = mean_bias / mean_bias(2);

        MSE = [ Bhat_raw(:,2)  Bhat_simp(:,2)  Bhat_run(:,2) Bhat_sur(:,2) Bhat_sinc(:,2)] ;
        MSE = (MSE - B(2)).^2 ./ (B(2))^2;
        mean_MSE = mean((MSE),1);

        Vtrue = V(1,1);

        if doGLS
            Vbias = [ gls_VB_raw(2,2) gls_VB_simp(2,2) gls_VB_run(2,2) gls_VB_sur(2,2) gls_VB_sinc(2,2) ]' ;
            Vbias = (Vbias - Vtrue) / Vtrue;
        else
            Vbias = [ ols_VB_raw(2,2) ols_VB_simp(2,2) ols_VB_run(2,2) ols_VB_sur(2,2) ols_VB_sinc(2,2) ]' ;
            Vbias = (Vbias - Vtrue) / Vtrue;
        end



        % Power calculation given mean, std.dev, and known effect size
        % we'll test the power of estimation of the second B param. in all
        % three types of differencing
        %     q = zeros(5,3);
        %     p = ones(5,3);
        %     Spower=zeros(5,3);
        q = zeros(5,2);
        p = ones(5,2);
        Spower=zeros(5,2);
        alpha = 0.05;
        alpha = 0.001;


        df = length(X) - length(B);

        % repeat for each differencing method
        for Dtype=1:5
            % repeat for each of the Bhats
            %         for c=1:3
            for c=1:2
                effect_size = abs(B(c));
                %Btmp = mean_Bhat(Dtype, c);
                sigma_B = std_Bhat(Dtype,c);

                % find which value of beta_hat corresponds to the
                % significance level alpha
                % in the null dustribution:
                tcrit =  spm_invTcdf(1-alpha, df);
                bcrit = tcrit * sigma_B;
                q(Dtype,c) = spm_Ncdf(bcrit, effect_size, sigma_B^2);

            end
            bcrit_SNR(Dtype,nl) = bcrit;
        end
        Spower = 1-q

        % and of course over several noise levels:
        %     Spower_SNR (:,nl) = Spower(:,3);
        %     eff_SNR (:,nl) = eff(:,3);
        %     bias_SNR (:,nl) = mean_bias(:,3);
        %     std_SNR (:,nl) = std_Bhat(:,3);
        Spower_SNR (:,nl) = Spower(:,2);
        eff_SNR (:,nl) = eff(:,2);
        bias_SNR (:,nl) = mean_bias;
        std_SNR (:,nl) = std_Bhat(:,2);

        % in the case of a monte carlo simulation, we'll average the power
        % and efficiency over iterations.  This means that we are not
        % keeping track of the noise level any more!
        Spower_MC (:,iter) = Spower(:,2);
        eff_MC (:,iter) = eff(:,2);

        % make a little table
        %     fprintf('\n Power  Efficiency   \n');
        %     table = [Spower(:,3)  eff(:,3)  ]
        %     table ./ repmat(table(1,:),5,1)

    end


    fprintf('\n Power  Efficiency  Variance bias mesn_MSE \n');
    table = [Spower(:,2)  eff(:,2)  Vbias mean_MSE']
    table ./ repmat(table(1,:),5,1)

end

mean_Spower_MC= mean(Spower_MC,2);
mean_eff_MC = mean(eff_MC,2);
std_Spower_MC= std(Spower_MC,0,2);
std_eff_MC = std(eff_MC,0,2);

if showPlots==1
    figure
    plot(1./nlevels, 100*Spower_SNR')
    legend('D1','D2','D3','D4','D5')
    xlabel('SNR'), ylabel('Power (%)')
    title('Power')
    axis tight, fatlines, dofontsize(16)

    figure
    plot(1./nlevels, eff_SNR);%'.*repmat(nlevels',1,4))
    %legend('D1','D2','D3','D4')
    xlabel('SNR'), ylabel('Efficiency')
    legend('D1','D2','D3','D4','D5')
    axis tight, fatlines, dofontsize(16)

%     figure
%     plot(1./nlevels, bias_SNR')
%     %legend('D1','D2','D3','D4')
%     xlabel('SNR'), ylabel('BIAS')
%     legend('D1','D2','D3','D4','D5')
%     axis tight, fatlines, dofontsize(16)
% 
%     figure
%     plot(1./nlevels, std_SNR);%'./repmat(nlevels',1,4))
%     %legend('D1','D2','D3','D4')
%     xlabel('SNR'), ylabel('Estimate Variance')
%     legend('D1','D2','D3','D4','D5')
%     axis tight, fatlines, dofontsize(16)

end

if doGLS==1
    str = sprintf('save %s_GLS.mat', ResultsFile);
else
    str = sprintf('save %s_OLS.mat', ResultsFile);
end
eval(str);

return

%some extra code to make the final figures.
load Block_sim_results_OLS.mat
set(0,'DefaultFigurePosition',[1 400 600 400])

subplot(121), hold on
plot(1./nlevels, 100*Spower_SNR')
subplot(122), hold on
plot(1./nlevels, eff_SNR')

load Block_sim_results_GLS.mat
subplot(121)
plot(1./nlevels, 100*Spower_SNR(1,:)', 'k')
subplot(122)
plot(1./nlevels, eff_SNR(1,:)', 'k')

subplot(121)
xlabel('SNR'), ylabel('Power (%)')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.', 'no sub + GLS')
legend('Location', 'SouthEast'), legend boxoff
title('Statistical Power')
axis ([0 1 0 110]), fatlines, dofontsize(14)
subplot(122)
xlabel('SNR'), ylabel('Efficiency')
title('Efficiency')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.', 'no sub + GLS','Location', 'SouthEast'), legend boxoff
fatlines, dofontsize(14)

print -dpng Block_meancenter

%%%% check the Monte Carlo results at Noise = 5
% order of columns is:  raw, simp, run, sur, sinc
% first the OLS
load MCsim_results_OLS

SPow = mean_Spower_MC;
Effs = mean_eff_MC;
std_SPow = std_Spower_MC;
std_Effs = std_eff_MC;

% then the GLS
load MCsim_results_GLS

SPow = [mean_Spower_MC(1); SPow]
Effs = [mean_eff_MC(1); Effs]

std_SPow = [std_Spower_MC(1); std_SPow]
std_Effs = [std_eff_MC(1); std_Effs]

figure
subplot(211), bar(SPow/SPow(1)), title('Relative Statistical Power')
hold on, errorbar([1:6], SPow/SPow(1), std_SPow/SPow(1),'o')
fatlines, dofontsize(16)
subplot(212), bar(Effs/Effs(1)), title('Relative Efficiency')
hold on, errorbar([1:6], Effs/Effs(1), std_Effs/Effs(1),'o')
fatlines, dofontsize(16)

print -dpng ERpower_eff_MC

%%%%%%%%%%%%%%%%%%%%%%%%
load Blocksim_results_OLS
% remember nlevels(6) == 5, NITER==1
SPow = Spower_SNR(:,6);
Effs = eff_SNR(:,6);

% then the GLS
load Blocksim_results_GLS

SPow = [Spower(1,6); SPow]
Effs = [eff_SNR(1,6); Effs]

figure
subplot(211), bar(SPow/SPow(1)), title('Relative Statistical Power')
hold on, errorbar([1:6], SPow/SPow(1), std_SPow/SPow(1),'o')
fatlines, dofontsize(16)
subplot(212), bar(Effs/Effs(1)), title('Relative Efficiency')
hold on, errorbar([1:6], Effs/Effs(1), std_Effs/Effs(1),'o')
fatlines, dofontsize(16)

print -dpng Blockedpower_eff_MC