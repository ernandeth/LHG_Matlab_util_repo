function B = biot(current, wire, loc)
%function B = biot(current, wire, loc)



dx=0.2;
MU = 1.26e-6;  % Permeability H/m

dwire = diff(wire,1);
B=0;

for r0=1:length(dwire) 
    B = B + ...
        ( cross(dwire(r0,:), (loc - wire(r0,:))) ...
        / norm(loc - wire(r0,:) )^3  );
end


B = B*MU*current/(4*pi);


return