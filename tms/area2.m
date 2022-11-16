%function result = area2()

%
% area2.m :

	name = input('Enter file name with points: ', 's');
	data = read_mat(name,4);

	tri = delaunay(data(:,1), data(:,3) )
pause
	trimesh(tri,data(:,1), data(:,3), data(:,2) )

%return
