function tmsplot2()


% Luis Hernandez
% Last Edit 10-8-97
% Input : four column ascii data file 
% Plots 3D intensity data over a surface.


	disp3d
	hold off
	hidden off
	count = 1;
	name = 'junk';

	while ~isempty(name)
		name = input('Select TMS File - Use ENTER to abort: ','s');
		if ~isempty(name),
			set(gca,'DataAspectRatio',[1 1 1]);
			displayfile(name,count);
			count = count +1;
			hold on;
		end;
	end;

return;

%%%%%%%%%%%%%%%%%%%%%%

function displayfile(name, color)
% put scaled data onthe screen

	
	fp = fopen(name);
	temp = fscanf(fp,'%f');
	sz = size(temp);
	fclose (fp);
	
	for i=1:sz(1)/4
		for j=1:4
			temp2(i,j) = temp((i-1)*4+j); 
		end
	end
	data_mat = temp2;
	
	sz = size(data_mat);
	
	y = data_mat(1:sz(1),1);
	z = data_mat(1:sz(1),2);
	x = data_mat(1:sz(1),3);
	c = data_mat(1:sz(1),4);
	
	xfid = zeros(100,1);
	yfid = zeros(100,1);
	zfid = zeros(100,1);
	j=1;
	
	% Make all topographical points have zero activation intensity.
	for i=1:size(c)
		if c(i) == -5 
			xfid(j) = x(i);
			yfid(j) = y(i);
			zfid(j) = z(i);
			j= j+1;
			c(i)=0;
		end;
	end;
	
	xfid = xfid(1:j-1);
	yfid = yfid(1:j-1);
	zfid = zfid(1:j-1);

	switch color
	case 1,
		s = 'b*';
	case 2,
		s = 'g*';
	case 3,
		s = 'm*';
	case 4,
		s = 'y*';
	otherwise,
		s = 'c*';
	end;
		plot3(x,y,z,s) ;

	
	hold on;

	plot3(xfid,yfid, zfid, 'r*');	

	% This is where the function 'junk' used to go 
return;	


%%%%
function j = junk()
	yesno = input('Is there an additional topographical data set?(y/n) :','s');
	if yesno =='y'

	% read in topographical data file
	% Data are concatenated to earlier topographical data

		name = input('Select Surface file (3 columns of data): ', 's');
		fp = fopen(name);
		temp = fscanf(fp,'%f');
		sz = size(temp);
		fclose (fp);
		
		for i=1:sz(1)/3
			for j=1:3
				temp2(i,j) = temp((i-1)*3+j); 
			end
		end
		data_mat = temp2;
		
		sz = size(data_mat);
		
		y2 = data_mat(1:sz(1),1);
		z2 = data_mat(1:sz(1),2);
		x2 = data_mat(1:sz(1),3);
		c2 = zeros(sz(1),1);
	
		x=cat(1,x,x2);
		y=cat(1,y,y2);
		z=cat(1,z,z2);
		c=cat(1,c,c2);
	
	end;
	
	% create mesh:
		
	xlin = linspace(min(x), max(x), 30);
	ylin = linspace(min(y), max(y), 30); 
	[X, Y] = meshgrid(xlin, ylin);
	Z = griddata(x,y,z,X,Y);
	C = griddata(x,y,c,X,Y);

	
	surf(X,Y,Z,C)
	colorbar
return;
