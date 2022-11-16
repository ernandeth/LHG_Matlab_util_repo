function deletecasfiles(cas)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: deletecasfiles.m 582 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% deletecasfiles.m - deletes list of files
% 
% INPUTS
% cas - cell array of strings, full filenames of files to be deleted
% 
% OUTPUTS
% [implicit] - files listed in cas will be deleted
% 
% EXAMPLE
% casFiles = {'c:\temp\trash.txt'; 'c:\temp\trash2.txt'};
% deletecasfiles(casFiles)

% Author - Krisanne Litinas
% $Id: deletecasfiles.m 582 2013-05-20 16:16:55Z klitinas $

cellfun(@delete,cas,'uni',false)
