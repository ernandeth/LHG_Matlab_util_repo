function [strDate,sdn] = getfiledate(strFile,fmt)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: getfiledate.m 590 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE:
% strDirDicom = getdicompath;
% casFilesDicom = getdicomfiles(strDirDicom);
% strFile = casFilesDicom{1};
% strDate = getfiledate(strFile)

if nargin == 1
    fmt = 29; % see datestr for fmt code
end
s = dir(strFile);
sdn = s.datenum;
strDate = datestr(sdn,fmt);
