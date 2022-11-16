function orbit(deg,stp,ax1);
%	Format:	 orbit(deg,stp);
%		rotates the camera deg amount every stp around the scene
%
%	Cengizhan
%	4.18.1997
%
%
[az el]=view;
rotvec=0:deg/stp:deg ;

for i=1:length(rotvec)
%	set(ax1,'CameraViewAngleMode','manual');
	view([az+rotvec(i) el])
	drawnow;
end;
