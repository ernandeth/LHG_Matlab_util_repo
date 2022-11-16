%The first section of this code is here for the purposes of exploring
%different values of AR+WN paramaters in order to match the power spectra
%that we found in real data.  Below that section is the section of code
%that emplots ASLnoise_eff_112905, which computes efficiency, power, etc
%for the different models.

addpath 'C:\Documents and Settings\Jeanette A Mumford\My Documents\Research\matlabcode'
addpath 'C:\Documents and Settings\Jeanette A Mumford\My Documents\Research\RAwork\ASL\Code'

%______________________________________AR+WN exploration_________________________
%trying to follow the AR+WN sturcture a little better
N = 258;        %Sample size
r=0.9;          %AR parameter
sigAR=0.02;     %AR variance
sigWN=2;      %WN variance

SNR=1;
V=AR1(r, N);    %make AR1 function is in the path...
V=(sigAR/(1-V(1,2)^2)).*V+sigWN.*eye(N);
%V=V/SNR^2;
%V=VD3;
Spec=abs(fft(full(V(:,floor(N/2))))).^2;

% Select first half - 
if rem(N,2),         % nfft odd
  select = (1:(N+1)/2+1)';  % don't include DC or Nyquist components
else
  select = (1:(N)/2+1)';
end

% Calculate the single-sided spectrum which includes the full power
Pwr = [2*Spec(select,:)];
Freq = (select-1)/max(select-1)*0.5;


% Pwr_V=Pwr;
% Freq_V=Freq;
% Pwr_D2=Pwr;
% Freq_D2=Freq;

%figure
plot(Freq,Pwr)
%axis([0, 0.4, 0, 550])

sqrt(V(1,1))
%Note:  I find that r=0.9, sigAR=0.1 and sigWN=10 is close to the resting
%voxel spectru

plot(Freq,Pwr_382, Freq, Pwr_6164, Freq, Pwr1_5411)
legend('SNR=0.382', 'SNR=0.6164', 'SNR-1.5411')

axis([0, max(Freq), 0, max(Pwr)]
text(0.4, 600, 'SNR=0.382')



%___________________________________Code to get power results______________________

clear

% Common stuff
doSpectra = 0;
showPlots = 0;
doMeanCenter = 0;
warning on
rho = 0.9;
VarAR = 0.02;
VarWN = 2;
TR = 1.4;
nyq = 1/(2*TR);
tlen = 258;
NITER = 1;
doNewMatrix = 1;  % generate a new design matrix or read it from a file?
%nlevels = 1./[0:0.2:2];
nlevels=1;

%VarAR/(1-rho^2)+VarWN
%with BOLD

% Blocked Design ASL
DesignType=-1;  % 1=blocked, 2=ER Fixed, 3=ER randomized
ResultsFile='Block_ASL';
%ASLnoise_eff_112905JM
ASLnoise_eff_060124JM






%%I'm going to make some efficiency plots!
%I know that changing the ratio of VarAR/VarWN and changing rho are the
%only things that make a difference in the bias formulae




DesignType=-3;
rho=0.9;
x=[0:0.1:1 1:0.5:12];
biasS=repmat(NaN, 6,length(x));
biasVarB=repmat(NaN, 6, length(x));
releff=repmat(NaN, 6, length(x));
for ind=1:length(x)
    
    VarWN=1;
    VarAR=x(ind);
    ASLnoise_eff_060124JM
    biasS(:,ind)=table(:,5);
    biasVarB(:,ind)=table(:,6);
    releff(:,ind)=eff_rel;
    
end
%since we redefined sigma_AR as VarAR/1-rho^2 I need to reflect that in the
%plot

xfig=x./(1-rho^2);
%figure
% 
% biasS5=biasS;
% biasVarB5=biasVarB;

subplot(2,1,1)
plot(xfig, 100*(biasS(1:5, :)'-1), 'linewidth', 2); 
xlabel('$\sigma^2_{AR}/\sigma^2_{WN}$', 'interpreter', 'latex','Fontsize', 20)

%ylabel('$\% \hspace{0.05 cm} \textrm{Bias}\hspace{0.05 cm} \hat{\sigma}^2_D$', 'interpreter','latex', 'Fontsize',20)
legend('None', 'Pairwise', 'Running', 'Surround', 'Sinc', 'location', 'east')
axis([0 25 -10 2])
set(gca, 'Fontsize', 15)

subplot(2,1,2)
plot(xfig, 100*(biasVarB(1:5, :)'-1),'linewidth', 2 );
xlabel('$\sigma^2_{AR}/\sigma^2_{WN}$', 'interpreter', 'latex','Fontsize', 20)

%ylabel('$\% \hspace{0.05 cm} \textrm{Bias}\hspace{0.05 cm} \widehat{\textrm{Var}}(c\hat{\beta})$', 'interpreter','latex', 'Fontsize', 20)
axis([0 25 -75 200])
set(gca, 'Fontsize',15)


%plot the relative efficiency

DesignType=-1;
rho=0.9;
x=[0:0.1:1 1:0.5:12];
releffBL=repmat(NaN, 6, length(x));
releffFER=repmat(NaN, 6, length(x));
releffRER=repmat(NaN, 6, length(x));

for ind=1:length(x)
    
    VarWN=1;
    VarAR=x(ind);
    
    DesignType=-1;
    ASLnoise_eff_060124JM
    releffBL(:,ind)=eff_rel;
    
    DesignType=-2;
    ASLnoise_eff_060124JM
    releffFER(:,ind)=eff_rel;
    
    DesignType=-3;
    ASLnoise_eff_060124JM
    releffRER(:,ind)=eff_rel;
    
end


xfig=x./(1-rho^2);
subplot(3,1,1)
plot(xfig, releffBL(1:5, :)', 'linewidth',2)
legend('None', 'Pairwise', 'Running', 'Surround', 'Sinc')
axis([0 25 0.75 1.05])
xlabel('$\sigma^2_{AR}/\sigma^2_{WN}$', 'interpreter', 'latex','Fontsize', 20)
ylabel('Relative Efficiency', 'Fontsize', 20)
set(gca, 'Fontsize',15)

subplot(3,1,2)
plot(xfig, releffFER(1:5, :)', 'linewidth',2)
axis([0 25 0.5 1.05])
xlabel('$\sigma^2_{AR}/\sigma^2_{WN}$', 'interpreter', 'latex','Fontsize', 20)
ylabel('Relative Efficiency', 'Fontsize', 20)
set(gca, 'Fontsize',15)

subplot(3,1,3)
plot(xfig, releffRER(1:5, :)', 'linewidth',2)
axis([0 25 0.5 1.05])
xlabel('$\sigma^2_{AR}/\sigma^2_{WN}$', 'interpreter', 'latex','Fontsize', 20)
ylabel('Relative Efficiency', 'Fontsize', 20)
set(gca, 'Fontsize',15)


%Now to vary rho

DesignType=-1;
VarAR=0.1;
VarWN=10;
x_rho=0.1:0.1:0.9;
biasS_rho=repmat(NaN, 6,length(x_rho));
biasVarB_rho=repmat(NaN, 6, length(x_rho));

for ind=1:length(x_rho)
  
    rho=x_rho(ind);
    ASLnoise_eff_060124JM
     biasS_rho(:,ind)=table(:,5);
    biasVarB_rho(:,ind)=table(:,6);
end

figure
subplot(2,1,1)
plot(x_rho, 100*(biasS_rho(1:5, :)'-1) )
xlabel('\rho')
ylabel('% Bias of \sigma^2-hat')
legend('None', 'Pairwise', 'Running', 'Surround', 'Sinc', 'location', 'east')
   
subplot(2,1,2)
plot(x_rho, 100*(biasVarB_rho(1:5, :)'-1) )
xlabel('\rho')
ylabel('% Bias of  Var(\beta)')     


%power plots
DesignType=-1; 
rho = 0.9;
VarAR = 0.02;
VarWN = 2;
OLSpowBLOCK=repmat(NaN,1,20);
GLSpowBLOCK=repmat(NaN,1,20);

for i=1:20
    
    nlevels=0.1*i
    ASLnoise_eff_060124JM
    OLSpowBLOCK(1,i)=table(2,1);
    GLSpowBLOCK(1,i)=table(6,1);
end

DesignType=-2; 
OLSpowFER=repmat(NaN,1,20);
GLSpowFER=repmat(NaN,1,20);

for i=1:20
    
    nlevels=0.1*i
    ASLnoise_eff_060124JM
    OLSpowFER(1,i)=table(2,1);
    GLSpowFER(1,i)=table(6,1);
    
end


DesignType=-3; 
OLSpowRER=repmat(NaN,100,20);
GLSpowRER=repmat(NaN,100,20);

for i=1:20
    
    nlevels=0.1*i
    
    for Rrep=1:100
    ASLnoise_eff_060124JM
    OLSpowRER(Rrep,i)=table(2,1);
    GLSpowRER(Rrep,i)=table(6,1);
    end
    
end


%can load the data from 'C:\Documents and Settings\Jeanette A Mumford\My Documents\Research\RAwork\ASL\Data'

load OLSpow_RER.txt
load GLSpow_RER.txt

OLSpow_mean=mean(OLSpow_RER);
GLSpow_mean=mean(GLSpow_RER);

OLSpow_sd=std(OLSpow_RER);
GLSpow_sd=std(GLSpow_RER);

perdiff_RER=mean((OLSpow_RER-GLSpow_RER)./GLSpow_RER);

perdiff_sd_RER=std((OLSpow_RER-GLSpow_RER)./GLSpow_RER);




del=0.1.*(1:20)./1.45;

%plot(del, 100*OLSpowBLOCK, 'y', del, 100*GLSpowBLOCK, 'b', del, 100*OLSpowFER, 'c', del, 100*GLSpowFER, 'm')


subplot(3,1,1)
plot(del, 100*OLSpowBLOCK, 'r', del, 100*GLSpowBLOCK, 'blue', 'linewidth', 2)
axis([0 1.3 0 100])
ylabel('Power', 'fontsize', 14)
xlabel('SNR', 'fontsize', 14)
tit=title('Block Design','fontsize', 14);
legend('OLS', 'GLS')
set(gca, 'Fontsize',15)

subplot(3,1,2)
plot(del, 100*OLSpowFER, 'r', del, 100*GLSpowFER, 'blue', 'linewidth', 2)
axis([0 1.3 0 100])
ylabel('Power', 'fontsize', 14)
xlabel('SNR', 'fontsize', 14)
tit=title('Event Related-Fixed ISI','fontsize', 14);
set(gca, 'Fontsize',15)

subplot(3,1,3)
plot(del, 100*OLSpow_mean, 'r',del, 100*GLSpow_mean, 'b',  'linewidth', 2)
hold on;
plot(del, 100*(OLSpow_mean+2*OLSpow_sd), 'r:', del, 100*(OLSpow_mean-2*OLSpow_sd), 'r:','linewidth',2)
plot(del, 100*(GLSpow_mean+2*GLSpow_sd),  'b:', del,100*(GLSpow_mean-2*GLSpow_sd), 'b:', 'linewidth', 2)
ylabel('Power', 'fontsize', 14)
xlabel('SNR', 'fontsize', 14)
tit=title('Event Related-Random ISI','fontsize', 14);
axis([0 1.3 0 100])
hold off;
set(gca, 'Fontsize',15)


plot(del, 100*(OLSpowBLOCK-GLSpowBLOCK)./(GLSpowBLOCK), 'c', del, 100*(OLSpowFER-GLSpowFER)./(GLSpowFER), 'm', del, 100*perdiff_RER, 'g',  'linewidth',2)
legend('Block design', 'ER-Fixed ISI', 'ER-Random ISI', 'location', 'east')
axis([0 1.3 -40 1])
ylabel('% Difference in Power')
xlabel('SNR')
hold on
plot(del, 100*(perdiff_RER+2*perdiff_sd_RER), 'g:',del, 100*(perdiff_RER-2*perdiff_sd_RER), 'g:','linewidth',2)
hold off





set(gca, 'Fontsize',15)
plot(del, OLSpowBLOCK, 'r', del, GLSpowBLOCK, 'blue')
plot(del, (OLSpowBLOCK-GLSpowBLOCK)./GLSpowBLOCK)

plot(del, OLSpowFER, 'r', del, GLSpowFER, 'blue')


plot(del, (OLSpowFER-GLSpowFER)./GLSpowFER)


%_______________________________________Plot for regressors _______________

isi = [1:10 20:30];
X = zeros(tlen,1);
X(round(isi))=1;
            X = conv(X,spm_hrf(TR));
            X = X(1:tlen);
            X = X/max(X);
            
Xbold=X;            
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

time=(1:40).*1.4;
          
subplot(311)
plot(time, Xbold, '-bo','MarkerFaceColor','b')
xlabel('time (seconds)')
tit=title('BOLD', 'fontsize', 14);
axis ([0 56 -0.2 1.2])


subplot(312)
plot(time, Xc, '-bo','MarkerFaceColor','b')
xlabel('time(seconds)')
tit=title('Baseline Perfusion','fontsize', 14);
axis ([0 56 -0.2 1.2])

subplot(313)
plot(time, X, '-bo','MarkerFaceColor','b')
xlabel('time(seconds)')
tit=title('Activation ASL signal','fontsize', 14);
axis ([0 56 -0.6 0.6])



