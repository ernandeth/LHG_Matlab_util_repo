useSingleAlpha=0;
useConfound = 1;
runEverything 
pairedTtest_spmJr
!cp -r /home/data/asl/groupResults /home/data/asl/grAlphas

useSingleAlpha=0;
useConfound = 0;
pairedTtest_spmJr
!cp -r /home/data/asl/groupResults /home/data/asl/grAlphasNoreg


%{
useSingleAlpha=1;
useConfound = 1;
runEverything
!cp -r /home/data/asl/groupResults /home/data/asl/grSingleAlpha

useSingleAlpha=1;
useConfound = 0;
pairedTtest_spmJr
!cp -r /home/data/asl/groupResults /home/data/asl/grSingleAlphaNoreg
%}

close all

figure; 
z1 = lightbox('/home/data/asl/grAlphas/groupResults/baseline/Zmap_0002',[-4 4],6);
figure; 
z2 = lightbox('/home/data/asl/grAlphasNoreg/groupResults/baseline/Zmap_0002',[-4 4],6);
figure; 
z3 = lightbox('/home/data/asl/grSingleAlpha/groupResults/baseline/Zmap_0002',[-4 4],6);
figure; 
z4 = lightbox('/home/data/asl/grSingleAlphaNoreg/groupResults/baseline/Zmap_0002',[-4 4],6);

close all

load /home/data/asl/grAlphas/groupResults/baseline/spmJr.mat
figure; 
h1 = hist(var_est(find(var_est)),[0:150]); 
e1 = sum(var_est(find(var_est)));

axis([0 140 0 3500])
load /home/data/asl/grAlphasNoreg/groupResults/baseline/spmJr.mat
figure; 
h2=hist(var_est(find(var_est)),[0:150]);
e2 = sum(var_est(find(var_est)));


axis([0 140 0 3500])
load /home/data/asl/grSingleAlpha/groupResults/baseline/spmJr.mat
figure; 
h3 = hist(var_est(find(var_est)),[0:150]);
e3 = sum(var_est(find(var_est)));

axis([0 140 0 3500])
load /home/data/asl/grSingleAlphaNoreg/groupResults/baseline/spmJr.mat
figure; 
h4 = hist(var_est(find(var_est)),[0:150]);
e4 = sum(var_est(find(var_est)));
axis([0 140 0 3500])
figure
plot([h1 ; h2; h3; h4]'); 
legend('alphas', 'alphas - no confound', 'single alpha', 'single alpha - no confound')
title('Histograms of Residuals')


close all

[all_subjs h]= read_img('/home/data/asl/grAlphas/groupResults/baseline/groupFlowEffects.img');
meanFlows = zeros(h.tdim,1);
for p=1:h.tdim
    tmp = all_subjs(p,:);
    tmp = tmp(isfinite(tmp));
    tmp = tmp(find(tmp));
    meanFlows(p)= mean(tmp);
end
f = meanFlows;

[all_subjs h]= read_img('/home/data/asl/grSingleAlpha/groupResults/baseline/groupFlowEffects.img');
meanFlows = zeros(h.tdim,1);
for p=1:h.tdim
    tmp = all_subjs(p,:);
    tmp = tmp(isfinite(tmp));
    tmp = tmp(find(tmp));
    meanFlows(p)= mean(tmp);
end
f = [f meanFlows];

subplot(211), plot(f);
legend('alphas', 'single alpha')
alpha = load('/home/data/asl/mfiles/alphaRegressor.txt')
%{
alpha = load('/home/data/asl/mfiles/alphas2b.txt');
a = mean(alpha,2); alpha = 1./a; 
alpha = alpha-mean(alpha);
alpha = alpha / max(alpha);
tmp = [alpha; alpha];
tmp(1:2:end) = alpha;
tmp(2:2:end) = alpha;
alpha = tmp;
save /home/data/asl/mfiles/alphaRegressor.txt alpha -ascii
%}
subplot(212); plot(alpha)
