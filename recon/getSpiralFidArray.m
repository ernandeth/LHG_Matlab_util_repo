%  Loads in the FID off of a spiral sequence and reorders the Fids.  In
%  other words, if slij corresponds to the ith slice and the jth 
%  interleave, Vnmrj natively stores fids in the following order:
%  [sl11 sl21 sl31...sln1 sl12 sl22...].  This function reshapes the data
%  such that interleaves corresponding to the same slice are ordered on the
%  2nd dimension and slices are ordered along the 3rd dimension.

function  [dta, procpar] = getSpiralFid(fileLocation)

% download the fid and parameter files.
[dta, procpar] = readvarian(fileLocation);

% Set up paramters
[n m] = size(dta);
sl = procpar.ns;   % Number of Slices
Nint = procpar.Nint;  % Number of Interleaves
end
