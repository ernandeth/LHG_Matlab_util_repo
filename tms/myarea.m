function result = area()

%
% area.m :
% Eliminates zero activation data points
% Computes center of of activation area
% Computes the surface area covered by a 3D set 
% of points by computing
% the areas of triangles formed by the different points.
%

global MAXX;

	total = 0.0;
	name = input('Enter file name with points: ', 's');
	data = read_mat(name,4);

%	% discard zero-activation data points
%	sz = size(data);
%	j=1;
%	for i=1 : sz(1)
%		if data(i,4) > 0,
%			temp(j,:) = data(i,:);
%			j=j+1;
%		end;
%	end;
%	data = temp
%
	
	% Compute center of activation	
	disp3d;
	center = ctr2(data);
	hold on;
	set(gca,'DataAspectRatio',[1 1 1]);
	plot3(data(:,3) ,data(:,1), data(:,2), 'g*');

	MAXX = max(data(:,4) );
	
	% Extract the indices of the verteces of the triangles
	% by Delaunay triangulation
	tri = delaunay(data(:,3), data(:,1) );
	sz = size(tri);
	ntrgl = sz(1);

	triangles = struct(...
	 		... %[1 ntrgl], ...
			'verts'	, zeros(3,3), ...
			'area '	, 0.0, ...
			'ctr'	, [0 0 0], ...
			'val'	, 0.0	 );

	% Create triangles 
	for i=1:ntrgl
		triangles(i).verts = tri_verts(data, tri, i);
		triangles(i).val = tri_val(data,tri,i); 
		triangles(i).area = tri_area_old(triangles(i).verts);
		triangles(i).ctr = tri_ctr(triangles(i).verts);

		% Draw triangles and 
		% add areas of all triangles		
		draw(triangles(i).verts, triangles(i).val);	
		total = total + triangles(i).area;
	end;
	result = total;


	Activation_Area = total
	Activation_Center = center	

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=draw(verts, c)

% Uses the patch function to draw all the delaunay triangles
% intensity is normalized to MAX
global MAXX

	v = zeros(3,3);
	% Extract the verteces from the matrix
	v(:,1) = verts(:,3);
	v(:,2) = verts(:,1);
	v(:,3) = verts(:,2);

	color=[c c c];
	if MAXX ~= 0
		color = color/MAXX; 
	else
		color = [0 0 1];
	end; 
	patch(v(:,1), v(:,2), v(:,3),color);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = tri_verts(data, tri, num)

% extracts the coordinates of each triangle from the triangle number
% and the original data

	result = zeros(3,3);
	for k=1:3
		result(k,:) = data(tri(num,k), 1:3);
	end;


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = tri_val(data, tri, num)

% Averages the stimulation values of the three verteces of each triangle

	intens = 0;
	for k=1: 3
		intens = intens + data(tri(num,k), 4)
	end;
	result = intens/3;

return

%%%%%%%%%%%%%%%%%%%%%%%%%
function result = tri_verts_old(data, num, N)

% tri_verts()
% Searches for the two closest points to each individual 
% point and returns a matrix containing the coordinates the 
% triangle formed by all three 

	p1 = data(num,:) ;
	result = zeros(3,3);
	close = 1000.0;
	closest = 1000.0;
	closenum = 1;
	closestnum = 1;

	for i = 1 : N
		if i ~= num ,
			p2 = data(i, :);
			dist = distance(p1,p2);

			if  (dist < closest),
				close = closest;
				closenum = closestnum;
				closest = dist;
				closestnum = i;		
			elseif (dist < close),
				close = dist;
				closenum = i;
			end;
		end
	end

	points = [num closenum closestnum];

	result = [	data(num,:);
			data(closenum, :);
			data(closestnum, :)];

%dummy = input('searched array for closest two points', 's');

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = tri_ctr(verts)

% tri_ctr()
% Calculates the center point of a triangle from the position of its
% verteces.
	result = zeros(1,3);

	result(1,1) = (verts(1,1) + verts(2,1) + verts(3,1) )/3;
	result(1,2) = (verts(1,2) + verts(2,2) + verts(3,2) )/3;
	result(1,3) = (verts(1,3) + verts(2,3) + verts(3,3) )/3;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result=distance(point1, point2)

% distance()
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = tri_area(v)
% this is an equation that I found in a matlab document.  i don't trust it 
% very much at all.

	result = 0.0;

	z1 = v(1,1); z2 = v(2,1); z3 = v(3,1);
	x1 = v(1,2); x2 = v(2,2); x3 = v(3,2);
	y1 = v(1,3); y2 = v(2,3); y3 = v(3,3);

	result = sqrt(((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))^2 + ...
			((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1))^2 + ...
			((x2-x1)*(y3-y1)-(y2-y1)*(y2-y1))^2)/2 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = tri_area_old(verteces)

% Computes the area of a triangle given its verteces in 3D space
% The triangle must be defined as a 3X3 matrix.  Each ROW is the XYZ
% coordinates of a vertex.

	% Extract the verteces from the matrix
	p1 = verteces(1,1:3);
	p2 = verteces(2,1:3);
	p3 = verteces(3,1:3);

	% Determine longest side of the triangle
	d = zeros(1,3);
	d(1) = distance(p1,p2);
	d(2) = distance(p2,p3);
	d(3) = distance(p1,p3);
		

	% Reorient triangle to fit general case
	% and translate triangle such that theta is at the origin

	a = d(1);
	b = d(3);
	c = d(2);
	x1 = p1;
	x2 = p3;	
	x3 = p2;

	% Translation:
	origin  = x1;
	x1 = x1 - origin;
	x2 = x2 - origin;
	x3 = x3 - origin;

	base = a;
	theta = acos( dot(x2,x3) / (b*a) );
	height = b*sin (theta);

	result = 0.5 * base * height;


return;
