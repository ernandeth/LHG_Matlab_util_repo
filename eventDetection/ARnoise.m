function [n]=ARnoise(T,rho,s2)
n=zeros(T,1);
for t=1:T,
    if(t==1)
        n(t)=randn*sqrt(s2)/sqrt(1-rho^2);
    else
        n(t)=rho*n(t-1)+sqrt(s2)*randn;
    end
end