function fsolvetest

H = 0.1*...
    [   -1     2       3
        4       -5     6
        7       8       -9];

C = -diag([1 2 3]);
Npoints = 500;

x = randn(3,Npoints);
u = randn(3,Npoints);

HC = [H C];
xu = [x;u];

dxdt = HC*xu;

dt = 0.1;

for n=2:Npoints
    %dxdt = H*x(:,n) + C*u(:,n);
    dxdt(:,n) = HC *xu(:,n);
    x(:,n) = x(:,n-1) + dxdt(:,n)*dt;
end

plot(x')


HC2 = dxdt * pinv(xu)

%HC0 = zeros(size(HC));
%HC2 = fsolve(@mysys, HC0)

return


function dxdt = mysys(HC)

residual = norm(dxdt - HC*xu);


return

