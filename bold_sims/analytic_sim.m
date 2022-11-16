t = linspace(0,8*pi);
f =  0.1.*sin(t);
plot(t,f),title('the function')
pause

ac = xcorr(f,'unbiased')';
ac = ac(size(ac)/2 : size(ac))

plot(ac),title('Autocorrelation')
pause

A = mean(f)

sn = sqrt(var(f));
se = 0;

rho = (ac - A^2) ./ (sn*sqrt(sn^2 + se^2));

plot(rho),title('Corr. COefficient vs. shift')
pause

for i=0:100
f1 = sin(t-0.2*i);
c = corrcoef(f,f1);
a(i+1)=c(2,1);
end

plot(a),title('Corr. COefficient vs. shift- the usual way')