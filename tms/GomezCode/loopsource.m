function [rs Jvec]=loopsource(center,normal,a,Npoints)
%[rs Jvec]=loopsource(center,a,theta,phi,deltaS)
%center-center of loop
%a=radius 
%theta- rotated about the x-z plane theta degrees 
%phi- a second rotation about the x' y' phi degrees
%rs- source locations
%Jvec-corresponding current sources(assuming clockwise current flowing)
if normal(1)==0 && normal(2)==0 && normal(3)==1
    A=eye(3);
else
A=rotmatrix(normal,[0 0 1]);
end
for i=1:Npoints
rs(i,:)=transpose(A*[a*cos(2*pi*(i-1)/Npoints); a*sin(2*pi*(i-1)/Npoints); 0]);
rs(i,1)=rs(i,1)+center(1);
rs(i,2)=rs(i,2)+center(2);
rs(i,3)=rs(i,3)+center(3);
end

for i=1:Npoints
Jvec(i,:)=(rs(i*(i~=Npoints)+1,:)-rs(i,:))/(2*pi*a);
end


