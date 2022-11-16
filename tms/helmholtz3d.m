function [E, Ex, Ey , Ez] = helmholtz(current, omega,  wire)
% function [E, Ex, Ey , Ez] = helmholtz(current, omega, wire)
%
% produces an E field using the Helmholtz equation over a cubic space 
% from -1 to 1 meters.
%
% wire - n x 3 matrix with coordinates of points along the wire
% the order of those points specifies the direction of the current (meters)
%
% omega - frequency of the current (Hz)
%
% current - amplitude of current in Amps a scalar
%

dx=0.05;

EPSILON = 8.854e-12; % permitivity e / V.m
MU = 1.26e-6;  % Permeability H/m

% compute the unit direction of the current vector
dwire = diff(wire,1);
vlengths = sqrt(sum(dwire.^2, 2));

[x,y,z] = meshgrid([-1:dx:1],[-1:dx:1],[-1:dx:1] ) ;

DIM= length(x);
x = reshape(x,DIM^3,1);
y = reshape(y,DIM^3,1);
z = reshape(z,DIM^3,1);

space = [x y z];

E = zeros(length(space),3);
k = omega* sqrt(MU*EPSILON);

% calculate each node's contribution to each point in space
for r=1:size(dwire,1)
    
    % distance from each point of wire to everywhere in space
    tmp = repmat(wire(r,:), length(space),1);
    dist = sqrt(sum( (space - tmp).^2 , 2) );
    
    G = exp(-i*k*dist) ./ dist; 

    fprintf('\rwire element: %d',r);
    
    % add the contribution of this node's E field to each point in space
    E = E + kron(G , dwire(r,:)/vlengths(r) );

end

E = -i*omega*MU*E;

Ex=reshape(E(:,1),DIM,DIM,DIM);
Ey=reshape(E(:,2),DIM,DIM,DIM);
Ez=reshape(E(:,3),DIM,DIM,DIM);

E = sqrt( Ex.^2 + Ey.^2 + Ez.^2);
return