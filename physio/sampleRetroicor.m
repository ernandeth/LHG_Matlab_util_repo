% sample script for physio correction
% you will need the files in this tar file and the ortho package, 
% which you should already have if you know waht's good for you
%
% step 1: read the physio  data from the scanner and change its format
% this generates a file called physio.dat

physdata = convertEXphysio('060224ds_phys_1',0.025)

% step 2: create a matrix with the physio data made up of 
% basis functions

PhysioMat = mkPhysioMat('physio.dat',0.025, 10, 36, 2);

% step 3:  estimate the parameters of the design matrix and remove iit from 
% the data .  This will generate a 4D image file called residuals.img.

rmReg('vol_blah_0001', PhysioMat);
  
% ------------------------------ %
% for ASL data, let's try to regress ou the CSF time course from signal.  
% Identify CSF voxels as the top 10% in the histogram

ASLPhysioMat = getASLPhysMat('vol_blah_0001');