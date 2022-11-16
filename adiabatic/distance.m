function result=distance(point1, point2)

% Computes the distance between two points in 3D space
% usage:	result=distance(point1, point2)
%
%	where point1 and point2 have to be row vectors.
%

	result = sqrt(...
		(point1(1) - point2(1))^2 + ...
		(point1(2) - point2(2))^2  + ...
		(point1(3) - point2(3))^2 );
return
