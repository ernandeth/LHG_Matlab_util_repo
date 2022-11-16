function  Y=rotx(X,angle,str);
%	Format : Y=rotx(X,angle,str);
% rotate along x axis by an angle alfa (in degrees clockwise when str='deg')
% in X every point set is a row
% cengizhan
% 4.5.96
%

if str=='deg',
	alfa=angle*pi/180;
else
	alfa=angle;
end;
Mx=[	1 	0	 			0; 
		0	cos(alfa)		-sin(alfa);
		0 	sin(alfa)		cos(alfa)]; 
		
Y=X*Mx;
