function removesubdirs(strDirTop,strPatternBasename)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: removesubdirs.m 616 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strDirTop = 'c:\temp\fmri_data\originals\s1369plas\pretwo\study_20091203';
% strPatternBasename = '^dcm$';
% removesubdirs(strDirTop,strPatternBasename);

% Get all subdirs
casDirsSub = getsubdirs(strDirTop);

% Get basenames
casBases = cellfun(@basename,casDirsSub,'uni',false);

% Get subset of subdirs where basename matches pattern
indRemove = findnonemptycells(regexp(casBases,strPatternBasename));

% Remove matching subdirs
for i = 1:length(indRemove)
    ind = indRemove(i);
    strDirSub = casDirsSub{ind};
    [SUCCESS,MESSAGE,MESSAGEID] = rmdir(strDirSub,'s');
    if ~SUCCESS
       warning('daly:common:fileio','could not remove "%s"',strDirSub) 
    end
end
