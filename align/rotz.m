function  Y=rotz(X,angle,str);
%
% rotate along z axis by an angle alfa (in degrees clockwise)
% in X every point set is a row
% cengizhan
% 4.5.96
%

if str=='deg',
	alfa=angle*pi/180;
else
	alfa=angle;
end;
Mz=[	cos(alfa) 	-sin(alfa)		0; 
		sin(alfa) 	cos(alfa)		0;
		0 		  	0				1]; 
		
Y=X*Mz;
