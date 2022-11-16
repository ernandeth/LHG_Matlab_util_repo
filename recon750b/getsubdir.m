function cas = getsubdir(strDir,strPat)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: getsubdir.m 599 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% use regexp to get subdirs that fit pattern under input dir
%
% EXAMPLE
% strDir = 'C:\temp\fmri_testing\c1316plas_control';
% strPat = '^study_\w+';
% cas = getsubdir(strDir,strPat)

[files,dirs]=spm_select('List',strDir,'.*');
casDirs = cellstr(dirs);
c = regexpi(casDirs, strPat, 'names');
iKeep = findnonemptycells(c);
cas = casDirs(iKeep);
