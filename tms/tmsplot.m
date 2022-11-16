function tmsplot()

% SUrface Plot of the skull from TMS positioning Pen data
% Luis Hernandez
% Last Edit 10-8-97
% Input : four column ascii data file 


name = input('Select  data filename to display: ','s');

tms3d

hold off
hidden off


fp = fopen(name);
temp = fscanf(fp,'%f');
sz = size(temp);
fclose(fp);

for i=1:sz(1)/3
	for j=1:3
		temp2(i,j) = temp((i-1)*3+j); 
	end
end
data_mat = temp2;

sz = size(data_mat);

x = data_mat(1:sz(1),1);
z = data_mat(1:sz(1),2);
y = data_mat(1:sz(1),3);
% c = data_mat(1:sz(1),4);



% simple Plotting:
plot3(x,y,z,'o');

set(gca,'DataAspectRatio',[1 1 1]);

xlabel('X');
ylabel('Y');
zlabel('Z');

axis manual;

pause;

hold on

% create mesh:
xlin = linspace(min(x), max(x), 30);
ylin = linspace(min(y), max(y), 30); 


[X, Y] = meshgrid(xlin, ylin);

Z = griddata(x,y,z,X,Y);
%C = griddata(x,y,c,X,Y);

%[X,Y,Z] = meshgrid(x,y,z);

surf(X,Y,Z)

