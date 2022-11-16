function casFiles = buildcasfiles(strDir,casFileList)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: buildcasfiles.m 578 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strDir = 'c:\path\of\interest';
% casFileList = {'fileOne.txt','fileTwo.dat','fileThree.log'};
% casFiles = buildcasfiles(strDir,casFileList)

casFiles = strcat(fixpath(strDir),filesep,casFileList)
