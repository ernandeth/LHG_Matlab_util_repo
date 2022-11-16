function M=getMy(rotx,roty,rotz)
%
%	gives the matrix R=XQyYZ for given angles
%	cengizhan 8.7.1996
%
cx=cos(rotx);
sx=sin(rotx);
cy=cos(roty);
sy=sin(roty);
cz=cos(rotz);
sz=sin(rotz);
Qy=[0 0 1;
	0 0 0;
	-1 0 0];
X=[	1 	0 	0; 
	0 	cx sx;
	0 	-sx cx];
Y=[	cy 	0 	sy; 
	0 	1 	0;
	-sy 0	cy];
Z=[	cz 	sz 	0; 
	-sz cz 	0;
	0 	0	1];
%M=X'*Qy'*Y'*Z';
M=X*Qy*Y*Z;
