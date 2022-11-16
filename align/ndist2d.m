%  A machine specific mex file, callback is :
%  ndist2D(X,width,height,ch,refno)
%		input argument ch=2 calculates the distancemap
%		input argument ch=4 saves the distancemapx (x is refno)
%		input argument ch=6 loads distancemapx to RAM
%		input argument ch=8 gives back the euclidean distance vector 
%			values of a 2 D array from RAM
%		input argument ch=9 gives back the euclidean distance vector 
			values of a 2 D array from DISK
%		input argument ch=10 clears up the distancemap from memory
%
% ADDITIONAL NOTES :
%	- Used for 2D space with sizes (width, height)
%	- refno is an integer will be the suffix after distancemap
%	- X is a n x 2 matlab matrix where each point is in a row
%	- NO negative coordinates or coordinates outside the plavne (will crash !)
%	
%		by Cengizhan Ozturk
%		9/1/1997
%
