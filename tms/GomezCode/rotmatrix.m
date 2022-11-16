function A=rotmatrix(new,old)
%rotmatrix(new,old)
%creates rotation transformation from old normal to new normal
u=cross(old,new)/norm(cross(old,new));
if new(1)==old(1) && new(2)==old(2) && new(3)==new(3)
 A=eye(3);
else
theta=acos(dot(new,old)/(norm(new)*norm(old)));
A(1,1)=(1-u(1)^2)*cos(theta)+u(1)^2;
A(1,2)=-u(3)*sin(theta)-u(1)*u(2)*cos(theta)+u(1)*u(2);
A(1,3)=u(2)*sin(theta)-u(1)*u(3)*cos(theta)+u(1)*u(3);
A(2,1)=u(3)*sin(theta)-u(1)*u(2)*cos(theta)+u(1)*u(2);
A(2,2)=(1-u(2)^2)*cos(theta)+u(2)^2;
A(2,3)=-u(1)*sin(theta)-u(2)*u(3)*cos(theta)+u(2)*u(3);
A(3,1)=-u(2)*sin(theta)-u(1)*u(3)*cos(theta)+u(1)*u(3);
A(3,2)=u(1)*sin(theta)-u(2)*u(3)*cos(theta)+u(2)*u(3);
A(3,3)=(1-u(3)^2)*cos(theta)+u(3)^2;
end