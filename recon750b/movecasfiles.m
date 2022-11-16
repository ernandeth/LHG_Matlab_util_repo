function movecasfiles(casFiles,strDirTo)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: movecasfiles.m 609 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% casFiles = {'c:\temp\trash.txt'; 'c:\temp\trash2.txt'};
% strDirTo = 'c:\temp\destdir';
% movecasfiles(casFiles,strDirTo)
casDirsTo = cellstr(repmat(strDirTo,length(casFiles),1));
cellfun(@movefile,casFiles,casDirsTo);
