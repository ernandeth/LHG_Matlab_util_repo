function [E, Ex, Ey , Ez]=Efield(dI_dt, wire, dims, FOV)
% function [Ex, Ey , Ez] = Efield(dI/dt, wire, dims)
%
% produces a E field using finite difference Time domain method
% over a cubic space from -1 to 1.
%
% wire - n x 3 matrix with coordinates of points along the wire
% the order of those points specifies the direction of the current
%
% dI/dt - a scalar representint the time derivative of the current through
% the wire (Amps /s)
%
% dims = a vector with the dimansions in voxels (i.e.- matrix size like 64 64 64))
%
% the equation is this  E = - mu / (4*pi) * dI/dt * integral( 1/|r -l| )dl
%
% field of view (m)

dr = [1 1 1]./dims;
dx=dr(1);dy=dr(2);dz=dr(3);  % meters

MU = 1.26e-6;  % Permeability H/m

% wire = [wire; wire(1, :)];  % close the loop
dwire = diff(wire,1);  
% find the midpoints of the segments
wire2 = wire(2:end,:); 
wire1 = wire(1:end-1,:);
wire = (wire1 + wire2)/2;

% normalize the direction vector
%ndwire = sqrt( dwire(:,1).^2  + dwire(:,2).^2 + dwire(:,3).^2);
%dwire = dwire ./ [ndwire ndwire ndwire];

mygrid = linspace(-FOV/2,FOV/2, dims(1));
[x,y,z] = meshgrid(mygrid,mygrid, mygrid ) ;

XDIM= length(x);
YDIM= length(y);
ZDIM= length(z);

x = reshape(x,XDIM*YDIM*ZDIM,1);
y = reshape(y,XDIM*YDIM*ZDIM,1);
z = reshape(z,XDIM*YDIM*ZDIM,1);

Ex = zeros(XDIM*YDIM*ZDIM,1);
Ey = zeros(XDIM*YDIM*ZDIM,1);
Ez = zeros(XDIM*YDIM*ZDIM,1);

space = [x y z];

E =zeros(length(space),3);
K = - dI_dt * MU /(4*pi);   % Amps/sec *H / m

% % calculate E for all points xyz
for r=1:size(space,1)
    xyz = space(r,:);
    xyz = repmat(xyz, size(wire,1),1);   
    R = sqrt(sum((xyz - wire).^2,2));
    R = repmat(R, 1, 3);
    integral = sum(dwire./R,1);
    E(r,:) = K * integral;
    %fprintf('\rvoxel: %d',r);        
end
% Efield_integral(Ex, Ey,Ez,wire,dwire,dI_dt,x,y,z)


% break down the Vecotr field into its components
% Ex=reshape(Ex,XDIM,YDIM,ZDIM);
% Ey=reshape(Ey,XDIM,YDIM,ZDIM);
% Ez=reshape(Ez,XDIM,YDIM,ZDIM);
 Ex=reshape(E(:,1),XDIM,YDIM,ZDIM);
 Ey=reshape(E(:,2),XDIM,YDIM,ZDIM);
 Ez=reshape(E(:,3),XDIM,YDIM,ZDIM);

E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
% This should be in volts / meter

return