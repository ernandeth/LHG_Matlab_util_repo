function newP=getnewP(R,T,P,C)
%
%	gives the p'= R*(P-C)+C+T for given parameters
%	cengizhan 8.7.1996
%
dum1(1,:)=P(1,:)-C(1);
dum1(2,:)=P(2,:)-C(2);
dum1(3,:)=P(3,:)-C(3);
dum2=R*(dum1);
CT=C+T;
newP(1,:)=dum2(1,:)+CT(1);
newP(2,:)=dum2(2,:)+CT(2);
newP(3,:)=dum2(3,:)+CT(3);
