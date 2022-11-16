function [x y z] = make_solenoid(Rad, len, Nturns, Npts)
% function [x y z] = make_solenoid(Rad, len, Nturns, Npts)
% 
dt = 0.1;
theta = linspace(0,2*pi*Nturns, Npts);
spiral=  Rad*exp(-i.*theta);



x = real(spiral) ;
y = imag(spiral);
z = linspace(-len/2, len/2, Npts);

x=x(:);
y=y(:);
z=z(:);

plot3(x,y,z);
axis square
fatlines
axis off
return
