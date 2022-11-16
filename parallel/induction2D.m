function volts = induction2D(dI_dt, source, target, showVectors)
% the units have to be in meteres and amps/sec.

% number of places where to sample the magnetic field
Nelements = 10;
x = linspace(target(1,1) , target(end,1), Nelements);
y = linspace(target(1,2) , target(end,2), Nelements);
z = zeros(size(x));
locations = [x; y ; z]';

MU = 1.26e-6;  % Permeability H/m
% this is the term that goes into Faraday's Law
K = -dI_dt * MU / (4*pi);
B = [0 0 0];

for r=1:Nelements

    % this is the vector from the source rung to the point in the coil
    r1 = source(1,:) - locations(r,:) ;
    % the current is along the z-axis
    dl1 = [0,0,1];
    D1 = sum((locations(r,:) - source(1,:)).^2)^0.5;
    B1 = cross(r1 , dl1) / (D1+eps)^2;


    % do the same thing for the second rung, whose current runs in the
    % opposite direction
    r2 = source(2,:) - locations(r,:) ;
    dl2 = [0,0,-1];
    D2 = sum((locations(r,:) - source(2,:)).^2)^0.5;
    B2 = cross(r2 , dl2) / (D2+eps)^2;

    if showVectors
        % option: visualize the B-field vectors induced by each rung:
        axis (0.3*[-1 1 -1 1])
        quiver(x(r), y(r), B1(1)/100, B1(2)/100,'r'); hold on; quiver(x(r), y(r), B2(1)/100, B2(2)/100,'g');  drawnow
    end
    
    % integrate the fields
    B = B + B1 + B2;

end

% this is the integral of the field on that surface (or line in this case)
B = B/Nelements;

% Compute the normal component of this vector relative to the
% target coil's plane?
targetplane = target(2,:) - target(1,:);
Bmag = norm(B);
Bproj = B * targetplane';
Bnormal = Bmag*sin(acos(Bproj/Bmag));

% use Faraday's law to scale it appopriately
volts = K * Bnormal;

return
