% function ASLnoise_eff_050801TN
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
% DesignType = 2;  % ASL:  1=blocked, 2=ER Fixed, 3=ER randomized
%                  % BOLD:-1=blocked,-2=ER Fixed,-3=ER randomized
% ResultsFile = 'Blocksim_results';
% doSpectra = 0;
% showPlots = 0;
% warning off
% rho = 0.3;
% VarR = 0;  % Ratio of White Noise variance to AR variance (0=AR only)
% TR = 1.4;
% nyq = 1/(2*TR);
% tlen = 258;
% NITER = 1;
% doNewMatrix = 1;  % generate a new design matrix or read it from a file?

% % determine a range of noise levels - SNRs..
% %nlevels = [0.5 1:15 20 30];
% %%%%%nlevels = [0.5:0.5: 4];
% Want SNR's from 0.1 thorugh 2
% SNR = 0.5./nlevels
%nlevels = 0.5./fliplr([0.1:0.1:0.3 0.5:0.25:2]);
%nlevels = 0.05./0.1;

%disp(['Working on : ' ResultsFile])


isASL = DesignType>0;

% Make the differencing matrices
D1 = eye(tlen);
D2 = zeros(tlen/2,tlen);
for count=1:tlen/2
    D2(count,count*2-1)=1;
    D2(count,count*2)=-1;
end
D3 = zeros(tlen-1,tlen);
for count=1:tlen-1
    D3(count,count)=(-1)^(count-1);
    D3(count,count+1)=(-1)^(count);
end
D4 = zeros(tlen-2,tlen);
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

nDtype = 6; % D1-D5 OLS plus D1 GLS
Spower_SNR = zeros(nDtype,length(nlevels));
eff_SNR = zeros(nDtype,length(nlevels));
bias_SNR = zeros(nDtype,length(nlevels));
bcrit_SNR = zeros(nDtype,length(nlevels));
Spower_MC = zeros(nDtype,length(nlevels),NITER);
eff_MC = zeros(nDtype,length(nlevels),NITER);
S2bias_MC = zeros(nDtype,NITER);
VBbias_MC = zeros(nDtype,NITER);

close all;

for nl=1:length(nlevels)
    fprintf('\rNoise level ...  %f ', nlevels(nl));
    noise_amp = nlevels(nl);
    simp = zeros(NITER, tlen/2);
    run = zeros(NITER,tlen-1);
    sinterp=zeros(NITER,tlen);
    orig=zeros(NITER,tlen);
    sincc=zeros(NITER,tlen);

    VB_raw = zeros(NITER,2);
    VB_simp = zeros(NITER,2);
    VB_run = zeros(NITER,2);
    VB_sur = zeros(NITER,2);
    VB_sinc = zeros(NITER,2);

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
            switch abs(DesignType)
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

            if isASL
                %%% Build ASL design mtx

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
            else
                %%% Build BOLD design mtx

                Xc = ones(tlen,0);
                Xb = ones(size(X));
            end

            % put it all together:
            % note that in the baseline case we use only the tag to prevent
            % too much colinearity with the activation vector
            X = [Xb Xc X];
        else
            % load /Users/hernan/matlab/flow/noise/ASLX.mat
            fprintf(' Reading design matrix ...')
            if isASL
                %load ASLrER.mat
                load /Users/hernan/matlab/flow/noise/ASLX3.mat
            else
                % load BOLDrER.mat
                load /Users/hernan/matlab/flow/noise/BOLDX.mat
            end
        end
        if isASL
            B = [100 -1 -0.5]';
            Bint = logical([0 0 1]);  % Beta's of interest
        else
            B = [100 1]';
            Bint = logical([0 1]);    % Beta's of interest
        end

        signal = X*B;

        % a better way to create AR1 noise?
        % note - now we vary SNR by changing the amplitude of the effect
        % Vo = spm_Q(rho,tlen);  % Correlation matrix
        Vo = WKfun('AR+WN',rho,VarR,tlen);
        % This is just a safety measure,to force homogeneous variance:
        Vo = WKfun('Cov2Cor',Vo);
        % Variance-covariance matrix
        V = Vo * noise_amp^2;
        % W = V^(-1/2) to make inversions easier down the line...
        W = WKfun('mkW',[],V);

        t=0:tlen-1;


        % method 1: no differencing
        D=D1;
        DX = D*X;
        Bi = Bint;
        zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = [];
        % recall that W = V^(-1/2)
        gls_VB_raw = pinv(W*DX)*W*D*V*D'*W'*(pinv(W*DX))';
        ols_VB_raw = pinv(  DX)*  D*V*D'   *(pinv(  DX))';
        DVD = WKfun('Cov2Cor',D*Vo*D');
        ols_S2bias_raw = (trace(DVD) - trace(DX'*DVD*DX*(DX'*DX)^(-1))) / ...
            (size(DX,1)-size(DX,2));
        ols_VBbias_raw = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ...
            ( Bi*ols_VB_raw*Bi' ) * ols_S2bias_raw ;


        % method 2: simple subtraction

        D=D2;
        DX = D*X;
        Bi = Bint;
        zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = [];
        % recall that W = (D*V*D')^(-1/2), and hence does not exist
        % for any differencing methods since D*V*D' is singular
        %%%gls_VB_simp = repmat(NaN,3,3);
        ols_VB_simp = pinv(DX)* D*V*D' *(pinv(DX))';
        DVD = WKfun('Cov2Cor',D*Vo*D');
        ols_S2bias_simp = (trace(DVD) - trace(DX'*DVD*DX*(DX'*DX)^(-1))) / ...
            (size(DX,1)-size(DX,2));
        ols_VBbias_simp = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ...
            ( Bi*ols_VB_simp*Bi' ) * ols_S2bias_simp;


        % method 3: running subtraction

        D=D3;
        DX = D*X;
        Bi = Bint;
        zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = [];
        % recall that W = (D*V*D')^(-1/2), and hence does not exist
        % for any differencing methods since D*V*D' is singular
        %%%gls_VB_run = repmat(NaN,3,3);
        ols_VB_run = pinv(DX)*D* V *D'*(pinv(DX))';
        DVD = WKfun('Cov2Cor',D*Vo*D');
        ols_S2bias_run = (trace(DVD) - trace(DX'*DVD*DX*(DX'*DX)^(-1))) / ...
            (size(DX,1)-size(DX,2));
        ols_VBbias_run = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ...
            ( Bi*ols_VB_run*Bi' ) * ols_S2bias_run;


        % method 4: surround subtraction

        D=D4;
        DX = D*X;
        Bi = Bint;
        zCol = find(all(DX==0)); DX(:,zCol) = []; Bi(zCol) = [];
        % recall that W = (D*V*D')^(-1/2), and hence does not exist
        % for any differencing methods since D*V*D' is singular
        %%%gls_VB_sur = repmat(NaN,3,3);
        ols_VB_sur = pinv(DX)*D* V *D'*(pinv(DX))';
        DVD = WKfun('Cov2Cor',D*Vo*D');
        ols_S2bias_sur = (trace(DVD) - trace(DX'*DVD*DX*(DX'*DX)^(-1))) / ...
            (size(DX,1)-size(DX,2));
        ols_VBbias_sur = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ...
            ( Bi*ols_VB_sur*Bi' ) * ols_S2bias_sur;


        % method 5: sinc subtraction

        D=D5;
        DX = D*X;
        Bi = Bint;
        zCol = find(all(abs(DX-0)<1e-6)); DX(:,zCol) = []; Bi(zCol) = [];
        % recall that W = (D*V*D')^(-1/2), and hence does not exist
        % for any differencing methods since D*V*D' is singular
        %%%gls_VB_sinc = repmat(NaN,3,3);
        ols_VB_sinc = pinv(DX)*D* V *D'*(pinv(DX))';
        DVD = WKfun('Cov2Cor',D*Vo*D');
        ols_S2bias_sinc = (trace(DVD) - trace(DX'*DVD*DX*(DX'*DX)^(-1))) / ...
            (size(DX,1)-size(DX,2));
        ols_VBbias_sinc = ( Bi*pinv(DX)*D*D'*(pinv(DX))'*Bi' *V(1,1)) / ...
            ( Bi*ols_VB_sinc*Bi' ) * ols_S2bias_sinc;


        var_Bhat = ...
            [ols_VB_raw(end,end);
            ols_VB_simp(end,end);
            ols_VB_run(end,end);
            ols_VB_sur(end,end);
            ols_VB_sinc(end,end);
            gls_VB_raw(end,end);
            ];

        S2bias = ...
            [ols_S2bias_raw;
            ols_S2bias_simp;
            ols_S2bias_run;
            ols_S2bias_sur;
            ols_S2bias_sinc;
            1
            ];

        VBbias = ...
            [ols_VBbias_raw;
            ols_VBbias_simp;
            ols_VBbias_run;
            ols_VBbias_sur;
            ols_VBbias_sinc;
            1
            ];

        eff = 1./var_Bhat;
        eff_rel = eff / eff(2);


        % Power calculation given mean, std.dev, and known effect size
        % we'll test the power of estimation of the second B param. in all
        % three types of differencing
        %     q = zeros(5,3);
        %     p = ones(5,3);
        %     Spower=zeros(5,3);
        q = repmat(NaN,nDtype,1);
        alpha = 0.05;
        alpha = 0.001;


        df = length(X) - length(B);

        % repeat for each differencing method
        effect_size = abs(B(end));
        for Dtype=1:nDtype

            sigma_B = var_Bhat(Dtype);

            % find which value of beta_hat corresponds to the
            % significance level alpha
            % in the null dustribution:
            tcrit =  spm_invTcdf(1-alpha, df);
            bcrit = tcrit * sqrt(sigma_B);
            q(Dtype) = spm_Ncdf(bcrit, effect_size, sigma_B);

            bcrit_SNR(Dtype,nl) = bcrit;
        end

        Spower = 1-q;

        % and of course over several noise levels:

        Spower_SNR (:,nl) = Spower;
        eff_SNR (:,nl) = eff;

        % in the case of a monte carlo simulation, we'll average the power
        % and efficiency over iterations.  This means that we are not
        % keeping track of the noise level any more!
        Spower_MC (:,nl,iter) = Spower;
        eff_MC (:,nl,iter) = eff;
        S2bias_MC (:,nl,iter) = S2bias;
        VBbias_MC (:,nl,iter) = VBbias;

        % make a little table
        %     fprintf('\n Power  Efficiency   \n');
        %     table = [Spower(:,3)  eff(:,3)  ]
        %     table ./ repmat(table(1,:),5,1)

    end


    table = [Spower eff Spower./Spower(end) eff./eff(end) S2bias VBbias];
    fprintf('\n    Power      Eff      RelPow    RelEff    S2bias    VBbias\n');
    disp(table)
end

if showPlots

    figure;  set(gcf,'Position',[1 1 250,450]);

    imagesc(X)
    colormap(gray)

    figure ,set(gcf,'Position',[1 1 450,450])
    plot(signal), title('signal'), drawnow
end

mean_Spower_MC= mean(Spower_MC,3);
mean_eff_MC = mean(eff_MC,3);
std_Spower_MC= std(Spower_MC,0,3);
std_eff_MC = std(eff_MC,0,3);

if showPlots==1
    figure
    plot(0.5./nlevels, 100*Spower_SNR')
    legend('D1','D2','D3','D4','D5')
    xlabel('SNR'), ylabel('Power (%)')
    title('Power')
    axis tight, fatlines, dofontsize(16)
    axis([0 0.5 0.5 1.1])

    figure
    plot(0.5./nlevels, eff_SNR);%'.*repmat(nlevels',1,4))
    %legend('D1','D2','D3','D4')
    xlabel('SNR'), ylabel('Efficiency')
    legend('D1','D2','D3','D4','D5')
    axis tight, fatlines, dofontsize(16)

    figure
    plot(0.5./nlevels, bias_SNR')
    %legend('D1','D2','D3','D4')
    xlabel('SNR'), ylabel('BIAS')
    legend('D1','D2','D3','D4','D5')
    axis tight, fatlines, dofontsize(16)

end

str = sprintf('save %s.mat', ResultsFile);

eval(str);

return

ResultsFile = 'Block_sim_results';
%ResultsFile = 'ERfixed_sim_results';
%ResultsFile = 'ERrand_sim_results';
%ResultsFile = 'ERrand_sim_resultsniter1';
RO = sprintf('%s_OLS',ResultsFile);
RG = sprintf('%s_GLS',ResultsFile);

load(RG)
GLSpower_SNR = Spower_SNR(1,:)';
GLSeff_SNR   = eff_SNR(1,:)';

%some extra code to make the final figures.
load(RO)
set(0,'DefaultFigurePosition',[1 400 1200 400])
figure

subplot(131)
plot(0.5./nlevels, 100*Spower_SNR'), hold on
subplot(132)
plot(0.5./nlevels, 100*Spower_SNR'./repmat(GLSpower_SNR,1,5)), hold on
subplot(133)
plot(0.5./nlevels, eff_SNR'), hold on

load(RG)
subplot(131)
plot(0.5./nlevels, 100*Spower_SNR(1,:)', 'k')
subplot(133)
plot(0.5./nlevels, eff_SNR(1,:)', 'k')

subplot(131)
xlabel('SNR'), ylabel('% Probability of Detection')
tit=title('Statistical Power');
axis ([0 2 0 100]), fatlines,  dofontsize(10)
hold off

subplot(132)
xlabel('SNR'), ylabel('Power as % of nosub+GLS')
tit=[tit title('Relative Power')];
fatlines, dofontsize(10)
hold off

subplot(133)
xlabel('SNR'), ylabel('1 / Var(c\beta)')
tit=[tit title('Efficiency')];
legend('no sub.','pairwise','running','surround','sinc', 'no sub GLS')
legend('Location', 'NorthWest'), legend boxoff
leg=legend;
fatlines, dofontsize(10)
hold off

set(tit,'fontsize',14)
set(leg,'FontSize',8)

set(gcf,'PaperPosition',[.5 1 [1200 400]/1200*7.5])
print -depsc2 ERrand.eps
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
plot(0.5/nlevels, 100*Spower_SNR' ./repmat(GLSpower_SNR,1,5))
xlabel('SNR'), ylabel('Power (%)')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.')
legend('Location', 'SouthEast'), legend boxoff
title('Statistical Power rel to GLS')
axis ([0 1 0 110]), fatlines, dofontsize(14)

subplot(122),
plot(0.5/nlevels, 100*eff_SNR' ./ repmat(GLSeff_SNR,1,5))
xlabel('SNR'), ylabel('Efficiency')
title('Efficiency rel to GLS')
legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.','Location', 'SouthEast'), legend boxoff
fatlines, dofontsize(14)

print -dpng Blocked



%%%% check the Monte Carlo results at Noise = 5
% order of columns is:  raw, simp, run, sur, sinc
% first the OLS
%load MCsim_results_OLS
load ERrand_sim_results_OLS

SPow = mean_Spower_MC;
Effs = mean_eff_MC;
std_SPow = std_Spower_MC;
std_Effs = std_eff_MC;

% then the GLS
%load MCsim_results_GLS
load ERrand_sim_results_GLS

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

print -dpng ERrand_bargraph

%%%%%%%%%%%%%%%%%%%%%%%%
load  VarianceCheck_ASL

SPow = Spower;
Effs =  eff;
std_SPow = std_Spower_MC;
std_Effs = std_eff_MC;


 
subplot(211), h=bar(SPow/SPow(1));set(h,'facecolor',[1 1 1]*2/3),
title('Relative Statistical Power')
hold on, errorbar([1:6], SPow/SPow(1), std_SPow/SPow(1),'linestyle','none','color','k');hold off
dofontsize(16);set(gca,'xlim',[0 7],'ylim',[0 1.2])
abline('h',1);
axis([0 7 0.8 1.1])
set(gca,'Xticklabel',{'none+OLS','pairwise','running','surround','sinc','none+GLS'})
subplot(212), h=bar(Effs/Effs(1));set(h,'facecolor',[1 1 1]*2/3),
title('Relative Efficiency')
hold on, errorbar([1:6], Effs/Effs(1), std_Effs/Effs(1),'linestyle','none','color','k');hold off
dofontsize(16); set(gca,'xlim',[0 7],'ylim',[0 1.2]); % fatlines,
abline('h',1,'color','k');
set(gca,'Xticklabel',{'none+OLS','pairwise','running','surround','sinc','none+GLS'})
axis([0 7 0.8 1.1])

print -depsc2 VarianceCheck_MC



