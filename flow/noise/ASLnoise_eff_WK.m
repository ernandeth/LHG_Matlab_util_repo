% From hernan@umich.edu Mon Jul 18 23:18:11 2005


% the batch script (comment out the top section of ASLnoise_eff)
clear

% Common stuff
doSpectra = 0;
showPlots = 0;
doMeanCenter = 0;
warning on
rho = 0.95;
VarR = 100;
TR = 1.4;
nyq = 1/(2*TR);
tlen = 258;
NITER = 1;
doNewMatrix = 1;  % generate a new design matrix or read it from a file?
nlevels = [0.5 1:15 20 30];



% Blocked Design ASL
DesignType=1;  % 1=blocked, 2=ER Fixed, 3=ER randomized
ResultsFile='Block_ASL';
ASLnoise_eff_050803TN

% ER fixed ASL
DesignType = 2;  % 1=blocked, 2=ER Fixed, 3=ER randomized
ResultsFile = 'ERFixed_ASL';
ASLnoise_eff_050803TN

% ER rand
DesignType=3;  % 1=blocked, 2=ER Fixed, 3=ER randomized
doNewMatrix = 0;  % generate a new design matrix or read it from a file?
NITER=1;
ResultsFile='ERrand_ASL';
ASLnoise_eff_050803TN
doNewMatrix = 1;

%%%
% doNewMatrix = 1;  % generate a new design matrix or read it from a file?
% NITER=100;
% ResultsFile='ERrand_sim_results2';
% ASLnoise_eff_050803TN


% Blocked Design BOLD
DesignType = -1;  % 1=blocked, 2=ER Fixed, 3=ER randomized
ResultsFile='Block_BOLD';
ASLnoise_eff_050803TN

% ER design BOLD
DesignType = -2;  % 1=blocked, 2=ER Fixed, 3=ER randomized
ResultsFile = 'ERFixed_BOLD';
ASLnoise_eff_050803TN

% ER rand BOLD
DesignType=-3;  % 1=blocked, 2=ER Fixed, 3=ER randomized
doNewMatrix = 0;  % generate a new design matrix or read it from a file?
NITER=1;
ResultsFile='ERrand_BOLD';
ASLnoise_eff_050803TN
doNewMatrix = 1;

prettyPlots2('Block_ASL')
prettyPlots2('ERFixed_ASL')
prettyPlots2('ERrand_ASL')
prettyPlots2('Block_BOLD')
prettyPlots2('ERFixed_BOLD')
prettyPlots2('ERrand_BOLD')


% checking the variance in the ER rand case
doNewMatrix = 1;
NITER = 100;
DesignType = 3;
ResultsFile = 'VarianceCheck_ASL'
nlevels = 5;
ASLnoise_eff_050803TN

