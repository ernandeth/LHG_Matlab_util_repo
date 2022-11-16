function xhat = deconv_tv_l1_fast(y,H,HtH,Hty,h1,h2,tolerance)
% Usage: xhat = deconv_tv_l1_fast(y,HtH,Hty,h1,h2,tol)
% estimates x from the cost function 
% J(x) = ||Hx-y||^2+h1*l1(x)+h2*TV(x)
% TV is the total variation penalty, l1 is the l1 penalty.
% tol is the tolerance
T=length(y);
rho=1e-8;
x0=y;
for i=1:1000
    delta_xt=diff(x0); 
    w1=h2/2./sqrt(delta_xt.^2+rho);
    W1=diag(w1);
    w=h1/2./sqrt(x0.^2+rho);
    W=diag(w);
    DtW1D=diag(w1(1:end))+diag([0; w1(1:end-1)])-diag(w1(1:end-1),1)-diag(w1(1:end-1),-1);
    DtW1D(T,T)=w1(end); DtW1D(T-1,T)=-w1(end); DtW1D(T,T-1)=-w1(end);
    A=HtH+W+DtW1D;
    x1=A\(Hty);
    res = y - H*x1;
    if( abs(x1-x0)<tolerance*T) break; end        
    x0=x1;
end
xhat=x1;