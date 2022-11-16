clear all
close all
clc

load DCM_data

x1 = X(:,1);
x2 = X(:,2);
[N,S]=nec_suf(x1,x2);
disp(['Probability of Necessity = ',num2str(N)])
disp(['Probability of Sufficiency = ',num2str(S)])
