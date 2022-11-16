function casDirs = getpatterndirs(strPattern,strDir)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: getpatterndirs.m 596 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% GETPATTERNDIRS use regular expression to find dirs & output as cas
%
% EXAMPLE
% strPattern = '^series.*mocoseries$';
% strDir = 'C:\temp\fmri_data\originals\c1367plas_two\study_20090922';
% casDirs = getpatterndirs(strPattern,strDir);

% Get dir(s)
cas = getsubdirs(strDir);

% Get matching dirs using strPattern as regexp
casBases = cellfun(@basename,cas,'uni',false);
indMatch = findnonemptycells(regexpi(casBases,strPattern));
casDirs = cas(indMatch);
