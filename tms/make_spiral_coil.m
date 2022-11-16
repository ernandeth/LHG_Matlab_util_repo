function [x y ] = make_spiral_coil(inrad, outrad, weight)
% function [x y ] = make_spiral_coil(inrad, outrad, weight)
% 
% this function generates a spiral coil that 
dt = 0.1;

Nturns = floor(weight);
Fturns = weight - Nturns;

r = linspace(outrad, inrad, 2*pi*abs(Nturns)/dt);
theta = linspace(0,2*pi*Nturns, length(r));
spiral=  (r.^3).*exp(-i.*theta);

if Fturns~=0
    r2 = outrad*sqrt(Fturns) * ones(1,2*pi/dt);
    theta2 = linspace(0,sign(weight)*(2*pi-0.1),length(r2));
    spiral2 = r2 .* exp(-i.*theta2);
    
    spiral = [spiral spiral2];
end

xshift = r.^3 - (3*outrad)*r ;

x = real(spiral) + xshift;
y = imag(spiral);
plot(x,y,'k');
axis square
fatlines
axis off
return
