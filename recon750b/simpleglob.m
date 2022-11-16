function cas = simpleglob(strPathPattern)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: simpleglob.m 617 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% GLOB   Filename expansion/search via wildcards.
%
% EXAMPLE
% strPathPattern = 'c:\temp\*\*ash*.txt';
% cas = simpleglob(strPathPattern)

% Look for python program on either unix or matlab path
% Unix first
[~,casTemp] = unix('which simpleglob.py');
casTemp = locparsesystemoutput(casTemp);
strPy = casTemp{1};

if ~isempty(strfind(strPy,'not found'))
    % if not on unix, try matlab path
    strPy = which('simpleglob.py');
    if isempty(strPy)
        error('fmrilab:utils:cmdnotfound','simpleglob.py not found, must be on either matlab and/or on unix path!');
    end
end

% Construct command and execute
strCmd = [strPy ' "' strPathPattern '"'];
[~,casTwo] = dos(strCmd);
cas = locparsesystemoutput(casTwo);

% ------------------------------------------
function casOut = locparsesystemoutput(casIn)
casOut = strsplit(10,casIn)';
indToss = findemptycells(casOut);
casOut(indToss) = [];
