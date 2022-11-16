function copycasfiles(casFiles,strDirTo)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: copycasfiles.m 581 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% 
% EXAMPLE
% casFiles = dirbs('s:\data\upper\clinical_measures\*hart*.xls');
% strDirTo = 'c:\temp\trash4kristen';
% copycasfiles(casFiles,strDirTo)

casDirsTo = cellstr(repmat(strDirTo,length(casFiles),1));
success = cellfun(@copyfile,casFiles,casDirsTo);

if ~all(success)
    fprintf('\nWarning: not all files copied successfully');
end