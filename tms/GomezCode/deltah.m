function dh=deltah(rs,Jvec,r0)
% dh=deltah(sourcelocation vector,Current vector, differentiallenght,observation point)
% returns the Biot-sabbart H=deltas/4*pi*JvecX|rs-r0|/|rs-r0|^3
%
R=[(r0(1)-rs(1)) (r0(2)-rs(2)) (r0(3)-rs(3))];
Rmag=sqrt(R(1)^2+R(2)^2+R(3)^2);
% if Rmag==0
% dh=[0 0 0];
% else
Rspec=R/Rmag^3;
dh=cross(Jvec,Rspec)/(4*pi);
% end