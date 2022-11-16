function bln = verifyfilecount(flist,strType,rng)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: verifyfilecount.m 623 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strPattern = '^a20.*bold.*\w+_(?<num>\d{1,2})3\.img$';
% rng = 11;
% strDir = pwd;
% strType = 'char'; % see getpatternfiles for strType(s)
% flist = getpatternfiles(strPattern,strDir,strType);
% bln = verifyfilecount(flist,strType,rng)

if isempty(flist)
    numFiles = 0;
else
    switch lower(strType)
        case 'char'
            numFiles = size(flist{1},1);
        case 'cas'
            numFiles = length(flist);
        otherwise
            error('wrong type for flist input')
    end
end

if ~isscalar(numFiles)
    error('numFiles not a scalar!?')
end

bln = ismember(numFiles,rng);

if ~bln
    fprintf('\nnumFiles = %d',numFiles)
    
end
