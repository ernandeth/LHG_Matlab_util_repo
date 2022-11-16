% leasqr example/test
global verbose
verbose=0;
t = [1:100]';
p=[1; 0.1];
data=leasqrfunc(t,p);
wt1=(1+0*t)./sqrt(data);  % = 1 /sqrt of variances of data
if (sscanf(version,'%f') < 4),
  rand('normal');
  data=data+0.05*rand(100,1)./wt1;   % 1/wt1 = sqrt of var = standard deviation
else
  data=data+0.05*randn(100,1)./wt1;
end;
options=[[0.01; 0.01] [.8; .8]];
[f1,p1,kvg1,iter1,corp1,covp1,covr1,stdresid1,Z1,r21]=...
  leasqr(t,data,[.8;.05],'leasqrfunc',.001,50,wt1,[1;1],...
  'leasqrdfdp',options);
