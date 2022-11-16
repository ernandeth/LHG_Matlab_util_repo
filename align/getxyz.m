function R=getXYZ(rotx,roty,rotz)
%
%	gives the matrix R=XYZ for given angles
%	cengizhan 8.7.1996
%
cx=cos(rotx);
sx=sin(rotx);
cy=cos(roty);
sy=sin(roty);
cz=cos(rotz);
sz=sin(rotz);
X=[	1 	0 	0; 
	0 	cx sx;
	0 	-sx cx];
Y=[	cy 	0 	sy; 
	0 	1 	0;
	-sy 0	cy];
Z=[	cz 	sz 	0; 
	-sz cz 	0;
	0 	0	1];
R=X*Y*Z;
