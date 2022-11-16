function  Y=roty(X,angle,str);
%
% rotate along y axis by an angle alfa (in degrees clockwise)
% in X every point set is a row
% cengizhan
% 4.5.96
%

if str=='deg',
	alfa=angle*pi/180;
else
	alfa=angle;
end;
My=[	cos(alfa) 	0	 	-sin(alfa); 
		0		 	1			0;
		sin(alfa) 	0		cos(alfa)]; 
		
Y=X*My;
