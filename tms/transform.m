function result = transform()

% Reads in TMS data file
% computes spatial transformation matrix from 
% two sets of fiducials using metrix.m
% Applies transformation matrix to data
% 
% USAGE: result = transform  (user will be prompted for file names)

	name = input('Select TMS data file to be transformed : ','s');
%	name = 'ch1cor.dat'

	data = read_mat(name)
	
	sz = size(data);

	xyz_data = data(1:sz(1), 1:3);
	xyz_data = [xyz_data  ones(sz(1),1)]';

	tms_data = data(1:sz(1), 4);

	A = metrix

	xyz_data = A*xyz_data;
	
	xyz_data = xyz_data';
	result = xyz_data(1:sz(1), 1:3);
	result = [result, tms_data]

	name = input('Type in file name for transformed output data: ', 's');
%	name = 'ch2.dat'

	write_mat(name, result);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = read_mat(name)
%  reads in a 4xN matrix from a text file:


	fp = fopen(name);
	temp = fscanf(fp,'%f');
	sz = size(temp);
	fclose (fp);

	for i=1:sz(1)/4
		for j=1:4
			mat(i,j) = temp((i-1)*4+j); 
		end
	end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = write_mat(name, data)
%  Writes 4xN matrix to a text file:

	fp = fopen(name, 'w');
	sz = size(data)
	
	for i=1:sz(1)
		for j=1:sz(2)
			fprintf(fp,'%f\t', data(i,j));
		end
		fprintf(fp,'\n');
	end

	fclose (fp);

return;
