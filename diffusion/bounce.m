function [x, y] = bounce(Nsteps, radius, scale)
%function [x, y] = bounce(Nsteps, radius, scale)
%
% this program returns the position in 2D of a particle moving inside a circle
% When the particle hits the wall, it bounces.
%
% Nsteps 	is the number of time steps
% radius 	is the radius of the circle (in whatever units)
% scale 	is the scale of the spatial dimensions of the jumps.  
% The steps are generated from 0 to 1.  Important that the radius and scale are in agreement
% 

%initialization of variables
x = zeros(size(Nsteps)); 
y = x; 
xx = x; yy=y; 
xvect = []; yvect = []; rvect = []; rxvect =[]; ryvect =[];
r=radius; 
verbose = 0;
warning off  % we ain't no wussies


for n = 2:Nsteps
     
    % gaussian size jumps.  The average size is the "scale"
    % which is the average jump computed from the eistein equation.
    % calculate the next position to jump to:
     x(n) = randn(1,1) * scale + x(n-1);
     y(n) = randn(1,1) * scale + y(n-1); 
     
    
     
     while (sqrt(x(n)^2 + y(n)^2)) > r 
         if (verbose >=1)
            fprintf( '\n outside circle ...%d\n',n)
         end
         
         %compute the slope of the line between the new point and the last point.
         m = (y(n) - y(n-1))/ (x(n)-x(n-1));
         b = y(n) - m*x(n);
       
         % calculate the point where that line intersects the circle
         A = 1+m^2;
         B = 2*m*b;
         C = b^2-r^2;
    
         X1test = (-B + sqrt(B^2 - 4*A*C))/(2*A);
         X2test = (-B - sqrt(B^2 - 4*A*C))/(2*A);
     
         Y1test = m*X1test + b;
         Y2test = m*X2test + b;
         
         r1test = sqrt(X1test^2 + Y1test^2);
         r2test = sqrt(X2test^2 + Y2test^2);
         
         
         if (verbose==2) 
            p = circle(r);
            hold on
            plot([X1test X2test] , [Y1test Y2test], 'y')
            drawnow
         end
            
         % decide which of the two solutions is the correct one
         % (closer to the end point)
         if (abs(x(n)- X1test)) < (abs((x(n) - X2test)));
             xx = X1test;
             yy = Y1test;
             
           else
             xx = X2test;
             yy = Y2test; 
             
             
         end 
         
         % calculate the vector of the movement and the 
         % tangent to the circle
         A = [x(n-1) y(n-1)];
         B = [xx yy];
         C = [0 0];
               
         A = A-B;
         B=B-B;
         C=C-B;
         % calculate the reflection angle
         theta = acos( dot(A,C) /  ( sqrt(A*A')*sqrt(C*C') )  );
         if isnan(theta)
             theta=0;
         end
         
         phi = atan(yy/xx);
         
         % calculate how much further to go after bounce
         dx = x(n) - xx;
         dy = y(n) - yy;
         dr = sqrt( dx^2 + dy^2);
         
         % calculate final position after bounce
         if xx < 0
            x2 = xx + dr*cos(theta + phi);
            y2 = yy + dr*sin(theta + phi);
        else   
            x2 = xx - dr*cos(theta + phi);
            y2 = yy - dr*sin(theta + phi);
        end
                
         
         if (verbose==2) 
            px = [x(n-1) ; xx ; x2];
            py = [y(n-1) ; yy ; y2]; 
            plot(px,py,'-b');
         end
         
         x(n) = x2;
         y(n) = y2;
         
      end
      
     
end    
if (verbose >= 1)
    hold off
    p = circle(r);
    hold on
    plot(x,y,'r')
    str = sprintf('Radius: %g', radius');
    title(str), drawnow
end      
return
