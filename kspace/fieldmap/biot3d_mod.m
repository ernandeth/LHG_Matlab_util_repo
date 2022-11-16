function [B, Bx, By , Bz]=biot3d_mod(current, wire)
% function [Bx, By , Bz] = biot3d(current, wire)
%
% produces a B field using the Biot-Savard Law over a cubic space 
% from -1 to 1.
%
% wire - n x 3 matrix with coordinates of points along the wire
% the order of those points specifies the direction of the current
%
% current - a scalar
%

dx=0.004;
MU = 1.26e-6;  % Permeability H/m

dwire = diff(wire,1);
wire = wire(2:end,:);

[x,y,z] = meshgrid([-0.128:dx:0.124],[-0.128:dx:0.124],[-0.128:dx:0.124] ) ; % changed to estimate Bfield only over the required FOV

DIM= length(x);
x = reshape(x,DIM^3,1);
y = reshape(y,DIM^3,1);
z = reshape(z,DIM^3,1);

space = [x y z];

B=zeros(length(space),3);




% for r=1:size(space,1)
%     
%     for r0=1:size(dwire,1) 
%         B(r,:) = B(r,:) + ...
%             ( cross(dwire(r0,:), (space(r,:)-wire(r0,:))) ...
%             / norm(space(r,:) - wire(r0,:) )^3  );
%     end
% 
% end

% faster way of doing it ...
% dwire = [dwire ; dwire(end,:)];
for r=1:size(space,1)
    
    sp = kron(space(r,:), ones(size(wire,1),1)  );
    
    num = cross(dwire ,(sp - wire), 2 );
    den =  sum( (sp - wire).^2 , 2 ).^1.5;
    den = [ den den den];
    
    B(r,:) = sum(num ./den,1);
%     fprintf('\rvoxel: %d',r);        
end

B = B*MU*current/(4*pi);

Bx=reshape(B(:,1),DIM,DIM,DIM);
By=reshape(B(:,2),DIM,DIM,DIM);
Bz=reshape(B(:,3),DIM,DIM,DIM);

return