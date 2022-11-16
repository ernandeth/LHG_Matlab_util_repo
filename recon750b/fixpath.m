function strPath=fixpath(str)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: fixpath.m 588 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% remove trailing filesep if it exists (for consistent path handling)
%
% strPath=fixpath(str);
%
% Inputs
% str - string for path with or without trailing filesep
%
% Outputs
% strPath - string same as str but without trailing filesep
%
% Example
% strSessionPath='/home/analysis2/robotics_data/s_rihafplas/RawData/20051202_Fri_pre/';
% fixpath(strSessionPath)

% $Author$ Hrovat
% $Id: fixpath.m 588 2013-05-20 16:16:55Z klitinas $

% Be nice about trailing filesep
strPath=str;
if strcmp(strPath(end),filesep)
    strPath(end)=[];
end
