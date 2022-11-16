function [t,p,est] = t_test(v1,v2)
% [t,p,est] = t_test(v1,v2)
% paired t-test
%

v1(isnan(v1)) = [];
v2(isnan(v2)) = [];

est = mean(v1) - mean(v2);

df = (length(v1) + length(v2) - 1);

sp = ((length(v1) - 1) * var(v1)) + ((length(v2) - 1) * var(v1)) ./ (df);
se = sp .* sqrt(1./length(v1) + (1./length(v2)));

t = est ./ se;



p = tdist(t,df);

return
