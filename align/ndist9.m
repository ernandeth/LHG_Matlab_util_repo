%  A machine specific mex file, callback is :
%  ndist9(X,width,height,depth,ch,refno)
%		input argument ch=2 calculates the distancemap
%		input argument ch=3 calculates and saves from disk (for low RAM)
%		input argument ch=4 saves the distancemapx (x is refno)
%		input argument ch=6 loads distancemapx to RAM
%		input argument ch=8 gives back the euclidean distance vector 
%			values of a 3 D array from RAM
%		input argument ch=9 gives back the euclidean distance vector 
			values of a 3 D array from DISK
%		input argument ch=10 clears up the distancemap from memory
%
% ADDITIONAL NOTES :
%	- Used for 3D space with sizes (width, height,depth)
%	- refno is an integer will be the suffix after distancemap
%	- X is a n x 3 matlab matrix where each point is in a row
%	- NO negative coordinates or coordinates outside the voxel space (will crash !)
%	
%		by Cengizhan Ozturk
%		9/1/1997
%
