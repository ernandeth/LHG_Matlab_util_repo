function [err] = LL_T1(T1_EST,tcurve, ti,alf, tau)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hh = size(ti);
numsam = hh(2);
magv = zeros(hh);

T1 = T1_EST(1);
M0 = T1_EST(2);

err = 0;

cosa = cos(alf);
sina = sin(alf);

beta = M0*(1-exp(-tau/T1))/((1-cosa*exp(-tau/T1))*sina);
T1s = tau/(tau/T1-log(cosa));
DR = (cosa*(1-(cosa*(exp(-tau/T1))^(numsam-1))))/((1+cosa*(cosa*(exp(-tau/T1)))^(numsam-1)))+1;   

for i = 1:hh(2)

         magv(i) = (beta*(1-DR*exp(-i*tau/T1s)));
       
end
err = sum((magv-tcurve).^2);

end

