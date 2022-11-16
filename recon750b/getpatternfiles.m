function casFiles = getpatternfiles(strPattern,strDir,strType)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: getpatternfiles.m 597 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% use regular expression to find files & output as cas
%
% EXAMPLE
% strPattern = '^a20.*bold.*\w+_(?<num>\d{2,3})\.img$';
% strDir = pwd;
% strType = 'cas'; % or 'char'
% casFiles = getpatternfiles(strPattern,strDir,strType);

% Get dir
strDir = [ fixpath(strDir) filesep ];

% Get matching files using strPattern as regexp
casFiles = {};
[files,dirs] = spm_select('List',strDir,strPattern);
if isempty(files), return, end

switch lower(strType)
    case 'cas'
        casFiles = strcat(strDir,cellstr(files));
    case 'char'
        casFiles = {strcat(strDir,files)};
end
