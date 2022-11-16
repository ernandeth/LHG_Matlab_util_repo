function result = tmsfit(data)

% Reads in TMS data file
% computes spatial transformation matrix from 
% two sets of fiducials using metrix.m
% Applies transformation matrix to data
% 
% USAGE: result = tmsfit  (data)

	sz = size(data);

	xyz_data = data(1:sz(1), 1:3);
	xyz_data = [xyz_data  ones(sz(1),1)]';

	tms_data = data(1:sz(1), 4);

	A = tmsfids;

	xyz_data = A*xyz_data;
	
	xyz_data = xyz_data';
	result = xyz_data(1:sz(1), 1:3);
	result = [result, tms_data]

return;

function A = tmsfids()

% Calculates the transformation matrix to make two coordinate
% systems match.  Assumes six fiducials in both the destination and 
% the source data
% USAGE : A = tmsfids
% Note:  assumes six fiducials in each set

	global to from;

   	[n p] = uigetfile('*.dat','Select file with destination fiducials');
   	to_name = strcat(p,n);
	[n p]= uigetfile('*.dat','Select file with source fiducials');
   	from_name = strcat(p,n);
   
	to = read_mat(to_name,4);
	from = read_mat(from_name,4);
	
	% Pad fiducial matrices with ones to make dimensions match 	
	to(:,4) = 1;
	from(:,4) = 1;

	to = to';
	from = from' ;

	% Initial guess for transformation matrix A
	A = ones(4,4);
	
	% Call least squares method to minimize difference between spaces
	A = leastsq('sse_fun',A);

	
	
return;


