p = load('MNI_points.txt');
plot3(p(:,1), p(:,2), p(:,3) ,'.')


surffit = fit([p(:,1), p(:,3)], p(:,2), 'poly23');
coeffs = coeffvalues(surffit);
sform = formula(surffit)

% desired surface y = f(x,z)
xmin = min(p(:,1));
xmax = - xmin;
zmin = min(p(:,3));
zmax = max(p(:,3))

[x z ] =  meshgrid(linspace(xmin, xmax, 50), linspace(zmin, zmax, 50)); 

y = surffit(x,z);
surf(x,y,z)
hold on
surf(1.1*x,1.1*y,1.1*z)
hold off