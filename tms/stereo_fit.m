function result = stereo_fit()

% 
% computes spatial transformation matrix from 
% two sets of fiducials 
% Applies transformation matrix to data
% 
% USAGE: result = stereofit

% scale the MRI fiducials to cm
[n p]=uigetfile('*.hdr', 'Select the image template');
h=read_hdr(strcat(p,n));
scale(1)=h.xsize; %(mm/pixel)
scale(2)=h.ysize;
scale(3)=h.zsize;

fids = read_mat('mri_fids.dat',3);
mm_fids(:,1) = fids(:,1) * scale(1);
mm_fids(:,2) = fids(:,2) * scale(2);
mm_fids(:,3) = fids(:,3) * scale(3);
write_mat(mm_fids,'mm_mri_fids.dat');

% compute the transformation matrix
disp ('Computing Affine Transformation ...')
A = Coregistration

data=[0 0 0];
xyz_data=[0 0 0]';


while size(data) ~= [0 0]
   
   data = input('Enter the pixel coordinates from the screen  [x y z]: ');
   
   if size(data) == [1 3]
      xyz_data = (data.*scale)'
      xyz_data(4) = 1	;			  % pad with a one
      result = A*xyz_data;
      disp('Coordinates in stereotaxic frame:  ')
      disp(result (1:3))
   end

end

return;

function A = Coregistration()

% Calculates the transformation matrix to make two coordinate
% systems match.  Assumes six fiducials in both the destination and 
% the source data
% USAGE : A = tmsfids
% Note:  assumes six fiducials in each set

	global to from;

   	to_name = 'stereo_fids.dat';
   	from_name = 'mm_mri_fids.dat';
   
	to = read_mat(to_name,3);
   from = read_mat(from_name,3);
   sz=size(to);
   
   to=[to ones(sz(1), 1)];
   from = [from ones(sz(1) , 1)];
		
	to = to';
	from = from' ;

	% Initial guess for transformation matrix A
	% A = eye(4,4);  
	A= [-0.9943    0.1264   -0.0069   86.0382
   -0.0255   -1.0440    0.2950   77.1127
    0.1467    0.1873    1.1717 -132.2902
         0         0         0    1.0000];

	% Call least squares method to minimize difference between spaces
	A = leastsq('sse_fun',A);

	
	
return;


