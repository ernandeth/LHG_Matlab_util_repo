function M=getRz2D(rotz)
%
%	gives the matrix R=QzX for a given angle
%	cengizhan 10.3.1996
%
cz=cos(rotz);
sz=sin(rotz);

M=[	cz 	sz; 
	-sz cz];

