% script to illustrate the properties of the different subtraction methods
% in ASL - turboCASL in event related experiments.
%
% uses simulated data with AR(1) noise
%
% 5 differencing matrices.
% calculates power, efficiency, MSE, bias of the beta estimate as a
% function of SNR
%

% clear
%
% % set up the options here:
% DesignType = 2;  % 1=blocked, 2=ER Fixed, 3=ER randomized
% ResultsFile = 'Blocksim_results';
% doSpectra = 0;
% showPlots = 0;
% doMeanCenter = 0;
% warning off
% doGLS = 1;  % other wise it's OLS
% rho = 0.3;
% TR = 1.4;
% nyq = 1/(2*TR);
% tlen = 258;
% NITER = 1;
% doNewMatrix = 1;  % generate a new design matrix or read it from a file?

% % determine a range of noise levels - SNRs..
% %nlevels = [0.5 1:15 20 30];
% %%%%%nlevels = [0.5:0.5: 4];
% nlevels = [0.5 20];

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

    Bhat_raw = zeros(NITER,2);
    Bhat_simp = zeros(NITER,2);
    Bhat_run = zeros(NITER,2);
    Bhat_sur = zeros(NITER,2);
    Bhat_sinc = zeros(NITER,2);

    VB_raw = zeros(NITER,2);
    VB_simp = zeros(NITER,2);
    VB_run = zeros(NITER,2);
    VB_sur = zeros(NITER,2);
    VB_sinc = zeros(NITER,2);

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
                    isi = [1:20:360 ]/1.4;
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

            % put in BOLD effect regressor here:

            % baseline signal
            Xb = ones(size(Xc));

            % put it all together:
            % note that in the baseline case we use only the tag to prevent
            % too much colinearity with the activation vector
            if (doMeanCenter)
                X = [Xc X];
                X = X-repmat(mean(X),tlen,1);
            else
                X = [Xb Xc X];
            end
        else
            %    load /Users/hernan/matlab/flow/noise/ASLX.mat
            fprintf(' Reading design matrix ...')
            load ASLX3.mat
        end
        if (doMeanCenter)
            B = [-1 -.5]';
        else
            B = [100 -1 -.5]';
        end
        signal = X*B;
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

        if doMeanCenter

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mean center the signal and noise cov V matrix here

            I = eye(size(V));
            J = ones(size(V))/length(signal);

            noise  = V^(1/2) *noise;
            raw    = noise + signal;

            raw    = (I-J)*raw;
            V      = (I-J)*V*(I-J)';    % Change V to reflect centering

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            noise = V^(1/2) *noise;
            raw   = noise + signal;
        end




        orig(iter,:) =  raw';
        t=0:tlen-1;


        % method 1: no differencing
        D=D1;
        DX = D*X; DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_raw (iter,:)= tmp(end-1:end)';
        if doGLS==1
            % recall that W = V^(-1/2)
            gls_VB_raw = pinv(W*DX)*D*W* V *W * D'*(pinv(W*DX))';
            VB_raw(iter,:) = (diag(gls_VB_raw(end-1:end,end-1:end)))';
        else
            ols_VB_raw = pinv(DX)*D* V *D'*(pinv(DX))';
            VB_raw(iter,:) = (diag(ols_VB_raw(end-1:end,end-1:end)))';
        end

        % method 2: simple subtraction

        xt = raw(1:2:end);
        xc = raw(2:2:end);

        asl = xt - xc;
        simp(iter,:) =  asl';

        D=D2;
        DX = D*X; DX(:,find(all(DX==0))) = [];

        tmp = pinv(DX)*D*raw;
        %        tmp = pinv(D*W*X)*D*W*raw;
        Bhat_simp (iter,:)= tmp(end-1:end)';

        if doGLS==1
            % recall that W = (D*V*D')^(-1/2), and hence does not exist
            % for any differencing methods since D*V*D' is singular
            gls_VB_simp = repmat(NaN,3,3);
            VB_simp(iter,:) = repmat(NaN,1,2);
        else
            ols_VB_simp = pinv(DX)*D* V *D'*(pinv(DX))';
            VB_simp(iter,:) = (diag(ols_VB_simp(end-1:end,end-1:end)))';
        end

        % method 3: running subtraction
        asl=zeros(tlen-1,1);
        for n=2:length(raw)
            asl(n-1) = (raw(n)-raw(n-1)) / 2 *(-1)^n;
        end
        run(iter,:)= asl';

        D=D3;
        DX = D*X; DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        Bhat_run(iter,:)= tmp(end-1:end)';

        if doGLS==1
            % recall that W = (D*V*D')^(-1/2), and hence does not exist
            % for any differencing methods since D*V*D' is singular
            gls_VB_run = repmat(NaN,3,3);
            VB_run(iter,:) = repmat(NaN,1,2);
        else
            ols_VB_run = pinv(DX)*D* V *D'*(pinv(DX))';
            VB_run(iter,:) = (diag(ols_VB_run(end-1:end,end-1:end)))';
        end

        % method 4: surround subtraction
        asl=zeros(tlen-2,1);
        for n=2:length(raw)-1
            asl(n-1) = (-raw(n+1)+ 2*raw(n) -raw(n-1)) / 4 * (-1)^n;
        end
        sur(iter,:)= asl';

        D=D4;
        DX = D*X; DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        %        tmp = pinv(D*W*X)*D*W*raw;
        Bhat_sur(iter,:)= tmp(end-1:end)';

        if doGLS==1
            % recall that W = (D*V*D')^(-1/2), and hence does not exist
            % for any differencing methods since D*V*D' is singular
            gls_VB_sur = repmat(NaN,3,3);
            VB_sur(iter,:) = repmat(NaN,1,2);
        else
            ols_VB_sur = pinv(DX)*D* V *D'*(pinv(DX))';
            VB_sur(iter,:) = (diag(ols_VB_sur(end-1:end,end-1:end)))';
        end

        % method 5: sinc subtraction
        D=D5;
        DX = D*X; DX(:,find(all(DX==0))) = [];
        tmp = pinv(DX)*D*raw;
        %        tmp = pinv(D*W*X)*D*W*raw;
        Bhat_sinc(iter,:)= tmp(end-1:end)';

        if doGLS==1
            % recall that W = (D*V*D')^(-1/2), and hence does not exist
            % for any differencing methods since D*V*D' is singular
            gls_VB_sinc = repmat(NaN,3,3);
            VB_sinc(iter,:) = repmat(NaN,1,2);
        else
            ols_VB_sinc = pinv(DX)*D* V *D'*(pinv(DX))';
            VB_sinc(iter,:) = (diag(ols_VB_sinc(end-1:end,end-1:end)))';
        end


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
                [diag(gls_VB_raw(end-1:end,end-1:end))';
                diag(gls_VB_simp(end-1:end,end-1:end))';
                diag(gls_VB_run(end-1:end,end-1:end))';
                diag(gls_VB_sur(end-1:end,end-1:end))';
                diag(gls_VB_sinc(end-1:end,end-1:end))';]
        else
            std_Bhat = ...
                [diag(ols_VB_raw(end-1:end,end-1:end))';
                diag(ols_VB_simp(end-1:end,end-1:end))';
                diag(ols_VB_run(end-1:end,end-1:end))';
                diag(ols_VB_sur(end-1:end,end-1:end))';
                diag(ols_VB_sinc(end-1:end,end-1:end))';
                ]

        end
        std_Bhat_rel = std_Bhat /std(Bhat_simp,0,1);

        eff = 1./std_Bhat;
        eff_rel = eff / eff(2);

        % bias and MSE on the second regressor (activation regressor)
        
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
                % we only want the 2nd and third regressors' estimates
                effect_size = abs(B(c+1));
                %Btmp = mean_Bhat(Dtype, c);
                sigma_B = std_Bhat(Dtype,c);

                % find which value of beta_hat corresponds to the
                % significance level alpha
                % in the null dustribution:
                tcrit =  spm_invTcdf(1-alpha, df);
                bcrit = tcrit * sigma_B;
                q(Dtype,c) = spm_Ncdf(bcrit, effect_size, sigma_B);

            end
            bcrit_SNR(Dtype,nl) = bcrit;
        end
        Spower = 1-q

        % and of course over several noise levels:
       
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

if showPlots

    figure;  set(gcf,'Position',[1 1 250,450]);

    imagesc(X)
    colormap(gray)

    figure ,set(gcf,'Position',[1 1 450,450])
    plot(signal), title('signal'), drawnow
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
    axis([0 0.5 0.5 1.1])

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

ResultsFile = 'Block_sim_results';
ResultsFile = 'ERfixed_sim_results';
ResultsFile = 'ERrand_sim_results';

% print -dpng Blocked

%%
%%  Now do it relatively
%%
%some extra code to make the final figures.
%load Block_sim_results_GLS.mat
load ERrand_sim_results_GLS
%load ERFixed_sim_results_GLS
GLSpower_SNR = Spower_SNR(1,:)';
GLSeff_SNR   = eff_SNR(1,:)';

%load Block_sim_results_OLS.mat
load ERrand_sim_results_OLS
%load ERFixed_sim_results_OLS
set(0,'DefaultFigurePosition',[1 400 600 400])

subplot(121),
plot(1./nlevels, 100*Spower_SNR' ./repmat(GLSpower_SNR,1,5))
xlabel('SNR'), ylabel('Power (%)')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.')
legend('Location', 'SouthEast'), legend boxoff
title('Statistical Power rel to GLS')
axis ([0 1 0 110]), fatlines, dofontsize(14)

subplot(122),
plot(1./nlevels, 100*eff_SNR' ./ repmat(GLSeff_SNR,1,5))
xlabel('SNR'), ylabel('Efficiency')
title('Efficiency rel to GLS')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.','Location', 'SouthEast'), legend boxoff
fatlines, dofontsize(14)

print -dpng Blocked



%%%% check the Monte Carlo results at Noise = 5
% order of columns is:  raw, simp, run, sur, sinc
% first the OLS
%load MCsim_results_OLS
load ERrand_sim_results_noise5_OLS

SPow = mean_Spower_MC;
Effs = mean_eff_MC;
std_SPow = std_Spower_MC;
std_Effs = std_eff_MC;

% then the GLS
%load MCsim_results_GLS
load ERrand_sim_results_noise5_GLS

SPow = [mean_Spower_MC(1); SPow]
Effs = [mean_eff_MC(1); Effs]

std_SPow = [std_Spower_MC(1); std_SPow]
std_Effs = [std_eff_MC(1); std_Effs]

std_SPow / SPow(1)
std_Effs / Effs(1)

figure
subplot(211), h=bar(SPow/SPow(1));set(h,'facecolor',[1 1 1]*2/3),
title('Relative Statistical Power')
hold on, errorbar([1:6], SPow/SPow(1), std_SPow/SPow(1),'linestyle','none','color','k');hold off
dofontsize(16);set(gca,'xlim',[0 7],'ylim',[0 1.2])
abline('h',1);
set(gca,'Xticklabel',{'none+GLS','none','pairwise','running','surround','sinc'})
subplot(212), h=bar(Effs/Effs(1));set(h,'facecolor',[1 1 1]*2/3),
title('Relative Efficiency')
hold on, errorbar([1:6], Effs/Effs(1), std_Effs/Effs(1),'linestyle','none','color','k');hold off
dofontsize(16); set(gca,'xlim',[0 7],'ylim',[0 1.2]); % fatlines, 
abline('h',1,'color','k');
set(gca,'Xticklabel',{'none+GLS','none','pairwise','running','surround','sinc'})

print -depsc2 ERrand_bargraph

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


