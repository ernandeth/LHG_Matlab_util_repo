function M=getMz2D(rotz)
%
%	gives the matrix R=QzX for a given angle
%	cengizhan 10.3.1996
%
cz=cos(rotz);
sz=sin(rotz);

Qz=[0 -1;
	1 0];
Z=[	cz 	sz; 
	-sz cz];

M=Qz*Z;
