function y = leasqrdfdp(x,f,p,dp,func)
%y = [0*x+1 x];
y=[exp(-p(2)*x) -p(1)*x.*exp(-p(2)*x)];
