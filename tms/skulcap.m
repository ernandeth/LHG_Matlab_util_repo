function skulcap(name)

% SUrface Plot of the skull from TMS positioning Pen data
% Luis Hernandez
% Last Edit 10-8-97
% Input : threee column ascii data file 
% Plots a mesh surface plot of the data after determining 
% relative position from reference point.

%hold off
%hidden off
tms3d

fp = fopen(name);
temp = fscanf(fp,'%f');
sz = size(temp);

for i=1:sz(1)/3
	for j=1:3
		temp2(i,j) = temp((i-1)*3+j); 
	end
end
data_mat = temp2;

sz = size(data_mat);

x = data_mat(1:sz(1),1);
y = data_mat(1:sz(1),2);
z = data_mat(1:sz(1),3);

for i=1 : sz(1)/2
	x1(i) = x(i*2) - x(i);
	y1(i) = y(i*2) - y(i);
	z1(i) = z(i*2) - z(i);
end

% simple Plotting:
%c = data_mat(1:sz(1),3);
%plot3(x1,y1,z1,'o');

set(gca,'DataAspectRatio',[1 1 1]);

% create mesh:
xlin = linspace(min(x1), max(x1), 50);
ylin = linspace(min(y1), max(y1), 50); 
[X, Y] = meshgrid(xlin, ylin);
Z = griddata(x1,y1,z1,X,Y);
C = ones(size(Z));

mesh(X,Y,Z,C)

