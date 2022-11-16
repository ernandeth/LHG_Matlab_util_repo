function A = metrix()

% Calculates the transformation matrix to make two coordinate
% systems match.  Assumes six fiducials in both the destination and 
% the source data
% USAGE : A = metrix
% Note:  assumes six fiducials in each set

	global to from;

	to_name = input('Type in destination fiducials file: ', 's');
	from_name = input('Type in source fiducials file: ','s');

%	to_name = 'ch1fids.dat'
%	from_name = 'ch2fids.dat'
	
	to = read_matrix(to_name)
	from = read_matrix(from_name)
	
	% Pad fiducial matrices with ones to make dimensions match 	
	to = [to ones(6,1)]';
	from = [from ones(6,1)]';

	% Initial guess for transformation matrix A
	A = ones(4,4);
	
	% Call least squares method to minimize difference between spaces
	A = leastsq('sse_fun',A);

	
	
return;


%%%%%%%%%%%%%%%%%%

function mat = read_matrix(name)
%  reads in a Nx3 matrix from a text file:


	fp = fopen(name);
	temp = fscanf(fp,'%f');
	sz = size(temp);
	fclose (fp);

	for i=1:sz(1)/3
		for j=1:3
			mat(i,j) = temp((i-1)*3+j); 
		end
	end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don Gage's notes:

%; from is a matrix containing the coordinates of the fiducials in the
%magnetic stimulator space

%from = [-8.385708568 1.051147237 10.43635383
%3.995189606 2.36526836 12.9643228
%-0.404430743 3.938382762 3.084894973
%0.205649109 -9.557679272 6.077234189]

%; bookkeeping details
%from = [from ones(4,1)]'

%; to is a matrix containing the coordinates of the fiducials in the MRI
%space

%to = [57 119 164
%211 99 129
%110 37 71
%127 142 11 ]

%; get to into same units as from
%to = to * 0.85938
%; bookkeeping details
%to = [to ones(4,1)]'


%; A is the matrix which will eventually hold the matrix for transforming
%the magnetic stimulation
%;   coordinates into MRI coordinates. The initialization with 1's is
%just to establish a starting point.
%;   I have not found it to be a problem in this application, but often
%finding the 'correct' answer and
%'   not a local minimum is very dependent on the starting guess (the A
%matrix initialization).

%A = [1 1 1 1
%1 1 1 1
%1 1 1 1
%1 1 1 1];

%; the following applies a nonlinear least squares approach for finding
%the elements of A. It uses
%;  the Levenberg-Marquardt method.  It is mimimizes the function sse_fun
%which is simply
%;  f = (to - A*from)^2. Hopefully it will find A such that to = A*from
%giving f = 0.

%A = leastsq('sse_fun',A);

