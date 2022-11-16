function [B, Bx, By , Bz]=biot3d(current, wire)
% function [Bx, By , Bz]=biot3d(current, wire)
%


dx=0.1;
MU = 1.26e-6;  % Permeability H/m

dwire = diff(wire,1);

[x,y,z] = meshgrid([-1:dx:1],[-1:dx:1],[-1:dx:1] ) ;

DIM= length(x);
x = reshape(x,DIM^3,1);
y = reshape(y,DIM^3,1);
z = reshape(z,DIM^3,1);

space = [x y z];

B=zeros(length(space),3);
B2=B;

wire2=wire(1:end-1, :);
tic
for r=1:size(space,1)
    
    for r0=1:size(dwire,1)
        num1= cross(dwire(r0,:), (space(r,:)- wire(r0,:)),2);
        den1=sum( (space(r,:) - wire(r0,:)).^2 )^1.5  ;
        B(r,:) = B(r,:) + num1/den1;
    end
 
end
toc
% other way of doing it ...
dwire = [dwire ; dwire(end,:)];
tic
for r=1:size(space,1)
    
    sp = kron(space(r,:), ones(size(wire,1),1)  );
    
    num = cross(dwire ,(sp - wire), 2 );
    den =  sum( (sp - wire).^2 , 2 ).^1.5;
    den = [ den den den];
    
    B2(r,:) = sum(num ./den,1);
        
end
toc


%B = B*MU*current/(4*pi);

Bx=reshape(B(:,1),DIM,DIM,DIM);
By=reshape(B(:,2),DIM,DIM,DIM);
Bz=reshape(B(:,3),DIM,DIM,DIM);

Bx2=reshape(B2(:,1),DIM,DIM,DIM);
By2=reshape(B2(:,2),DIM,DIM,DIM);
Bz2=reshape(B2(:,3),DIM,DIM,DIM);

return