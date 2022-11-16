
function result = rotation(data)
% result = rotation(data)
%
% Creates rotation matrix from rotation angles and applies it
% to the position coordinates.  
% The function takes in a six column matrix where the first
% three coordinates are the cartesian 3D coordinates and the 
% last three coordinates are the pitch roll and yaw of the 
% reference point.
% the function retruns a four column matrix (xyz + a column of zeros)
%
% Las Edit 1-5-98
% Luis Hernandez

	sz = size(data);
	result = zeros(sz(1),4);

	for i=1: sz(1),
		a = data(i,4) * 2*pi/360;
		e = data(i,5) * 2*pi/360;
		r = data(i,6) * 2*pi/360;

		rot = zeros(3,3);

		rot(1,1) = cos(e)*cos(a);
		rot(1,2) = cos(e)*sin(a);
		rot(1,3) = -sin(e); 
		rot(2,1) = -cos(r)*sin(a) + sin(r)*sin(e)*cos(a);
		rot(2,2) = cos(r)*cos(a) + sin(r)*sin(e)*sin(a);
		rot(2,3) = sin(r)*cos(e);
		rot(3,1) = sin(r)*sin(a) + cos(r)*sin(e)*cos(a);
		rot(3,2) = -sin(r)*cos(a) + cos(r)*sin(e)*sin(a);
		rot(3,3) = cos(r)*cos(e);

		result(i,1:3) = (rot * data(i,1:3)' )';
		rot;
	end

return