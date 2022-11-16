
alpha = load('/home/data/asl/mfiles/alphas2b.txt');
a = mean(alpha,2); 
alpha = 1./a;
alpha = alpha-mean(alpha);
alpha = alpha / max(alpha);

tmp = [alpha; alpha];
tmp(1:2:end) = alpha;
tmp(2:2:end) = alpha;
alpha = tmp;
save /home/data/asl/mfiles/alphaRegressor.txt alpha -ascii

RT1 = load('/home/data/asl/mfiles/RTime_1back.txt');
RT4 = load('/home/data/asl/mfiles/RTime_4back.txt');
P1 = load('/home/data/asl/mfiles/perform_1back.txt');
P4 = load('/home/data/asl/mfiles/perform_4back.txt');

tmp = RT4;

pre = tmp(:,1:2);
pre = pre'; pre = pre(:);
post = tmp(:,3:4);
post = post'; post = post(:);
tmp  = [pre; post];

save /home/data/asl/mfiles/RT4.txt tmp -ascii

%%
close all
figure;
subplot(211)
z1=lightbox('/home/data/asl/groupResultsAlphas/baseline/Zmap_0002',[-5 5],4);
load /home/data/asl/groupResultsAlphas/baseline/spmJr.mat
subplot(212)
hist(var_est(find(var_est)),100)

figure;
subplot(211)
z2=lightbox('/home/data/asl/groupResultsSingleAlpha/baseline/Zmap_0002',[-5 5],4);
load /home/data/asl/groupResultsSingleAlpha/baseline/spmJr.mat
subplot(212)
hist(var_est(find(var_est)),100)

figure;
subplot(211)
z3=lightbox('/home/data/asl/groupResultsSingleAlphaNoR/baseline/Zmap_0002',[-5 5],4);
load /home/data/asl/groupResultsSingleAlphaNoR/baseline/spmJr.mat
subplot(212)
hist(var_est(find(var_est)),100)

