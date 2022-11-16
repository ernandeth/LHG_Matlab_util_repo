function touchcasfiles(strDate,casFiles)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: touchcasfiles.m 622 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strDate = '1984-04-30';
% casFiles = {'C:\data\fmri\fromopticalmedia\temp4pre\s1374plas\13740000\dcm\68911178';'C:\data\fmri\fromopticalmedia\temp4pre\s1374plas\13740000\dcm\68911194'};
% touchcasfiles(strDate,casFiles)

casDates = cellstr(repmat(strDate,length(casFiles),1));
cellfun(@touch,casDates,casFiles);
