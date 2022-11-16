function result = tms_filt()
echo on
% TMS_FILT
% Reads in TMS data file
% computes spatial translation by subtracting coordinates from reference coordinates
% Discards the Separation Markers based on min and max thresholds
% Rotates data according to pitch, yaw, roll of reference marker
% Writes data to user selected location
% 
% USAGE: result = tms_filt  (user will be prompted for file names)
echo off

	name = input('Enter TMS RAW position data file name : ','s');
	output_name = input('Output data file name : ','s');

	data = read_mat(name, 6);
	
	input('translating data:...','s');

	data = translation(data);

maxi = 30;
mini = -30;

	input('Removing out of range markers...','s');
	data = remove(data, mini, maxi);


	%input('Selecting every fourth point','s');
	%data = reduce(data,4);

	input('Rotating data....','s');
	data = rotation(data);

	write_mat(output_name, data);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = reduce(data, ntrvl)
% Take every nth data point and throw out the rest.

	sz = size(data);
	result = zeros(sz(1)/4, sz(2) );
	for i=4 :4: sz(1)
		result(i/4, :) = data(i,:);
	end;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = rotation(data)
% create rotation matrix from rotation angles:

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
	end,

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = remove(data, mini, maxi)

	sz = size(data);
	rows = sz(1);
	cols = sz(2);
	j=1;

	for i=1:rows,
		if (data(i,1) < maxi)  & ( data(i,1) > mini) 
			result(j,:) = data(i,:);
			j=j+1;
		end;
	end;


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = translation(data)

% Performs the translation by subtractions cordinates from reference

	sz = size(data);
	rows = sz(1);
	cols = sz(2);
	result = zeros(rows/2,cols);
	for i= 1:rows/2,
		result(i,1:3) = data(i*2-1, 1:3) - data(i*2 , 1:3);
		result(i,4:6) = data(i*2,4:6);
	end;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = read_mat(name,cols)
%  reads in a 4xN matrix from a text file:


	fp = fopen(name);
	temp = fscanf(fp,'%f');
	sz = size(temp);
	fclose (fp);

	for i=1:sz(1)/cols
		for j=1:cols
			mat(i,j) = temp((i-1)*cols+j); 
		end
	end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = write_mat(name, data)
%  Writes NxN matrix to a text file:

	fp = fopen(name, 'w');
	sz = size(data);
	
	for i=1:sz(1)
		for j=1:sz(2)
			fprintf(fp,'%f\t', data(i,j));
		end
		fprintf(fp,'\n');
	end

	fclose (fp);
	input('Writing Output file...','s');
return;
