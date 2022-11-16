function y = leasqrfunc(x,p)
%sprintf('called leasqrfunc(x,[%e %e]\n', p(1),p(2))
%y = p(1)+p(2)*x;
y=p(1)*exp(-p(2)*x);
