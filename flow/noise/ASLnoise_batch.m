clear

doSpectra = 0;
showPlots = 0;
doMeanCenter = 0;
warning off
doGLS = 0;  % other wise it's OLS
rho = 0.9;
TR = 1.4;
nyq = 1/(2*TR);
tlen = 258;
NITER = 1;
doNewMatrix = 1;  % generate a new design matrix or read it from a file?
nlevels = [0.5 1:15 20 30];

% Blocked Design
DesignType=1;
ResultsFile='Block_sim_results';
doGLS=1;
ASLnoise_eff_050728TN
doGLS=0;
ASLnoise_eff_050728TN
prettyPlots(ResultsFile)

%  Fixed ER
DesignType = 2;  % 1=blocked, 2=ER Fixed, 3=ER randomized  
ResultsFile = 'ERFixed_sim_results';
doGLS=0;
ASLnoise_eff_050728TN
doGLS=1;
ASLnoise_eff_050728TN
prettyPlots(ResultsFile)

% rand ER Design (always do the same design)
doNewMatrix = 0;
DesignType=3;
ResultsFile='ERrand_sim_results';
doGLS=0;
ASLnoise_eff_050728TN
doGLS=1;
ASLnoise_eff_050728TN
prettyPlots(ResultsFile)

%{
% rand ER Design noise = 5, iterate to look at variance
doNewMatrix = 1;
DesignType=3;
ResultsFile='ERrand_sim_results_noise5';
NITER=100;
nlevels=[5];
ASLnoise_eff_050728TN
doGLS=0;
ASLnoise_eff_050728TN
prettyPlots(ResultsFile)
%}

% extract efficiencies table
E = zeros(3,6);
load Block_sim_results_GLS.mat               
E(1,1) = eff_SNR(1,6);
load Block_sim_results_OLS.mat
E(1,2:end) = eff_SNR(:,6)'


load ERFixed_sim_results_GLS.mat       
E(2,1) = eff_SNR(1,6);
load ERFixed_sim_results_OLS.mat  
E(2,2:end) = eff_SNR(:,6)'

load ERrand_sim_results_GLS.mat          
E(3,1) = eff_SNR(1,6);
load ERrand_sim_results_OLS.mat  
E(3,2:end) = eff_SNR(:,6)'

Erel = E./repmat(E(:,1),1,6)
