function strParentDir = pdir(strChildDir)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: pdir.m 610 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% return string that is parent dir of either pwd when no arg or of
% child input dir arg
%
% strParentDir = pdir(strChildDir);
%
% INPUTS:
% strChildDir (optional) - string for child dir (pwd if no input arg)
%
% OUTPUTS:
% strParentDir - string for parent dir
%
% EXAMPLE
% fprintf('\n\nThe parent dir of...\n%s\nis...\n%s\n\n',pwd,pdir)

% author: Ken Hrovat
% $Id: pdir.m 610 2013-05-20 16:16:56Z klitinas $

if nargin == 0
    strChildDir = pwd;
end
if ~exist(strChildDir,'dir')
    error('daly:common:fileNotFound','directory "%s" does not exist',strChildDir)
end
strParentDir = fileparts(strChildDir);
