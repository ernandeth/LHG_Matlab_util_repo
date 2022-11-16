function touch(strDate,strFile)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: touch.m 621 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strDate = '1984-04-30';
% strFile = 'C:\data\fmri\fromopticalmedia\temp4pre\s1374plas\13740000\dcm\68911178';
% touch(strDate,strFile)

%% Like this !c:\cygwin\bin\touch.exe -d "1984-04-30" C:\data\fmri\fromopticalmedia\temp4pre\s1374plas\13740000\dcm\68911178
strCmdTouch = ['c:\cygwin\bin\touch.exe -d "' strDate '" '];

%% Anonymize dicoms
[status,strResult] = system([strCmdTouch strFile]);
if max(size(strResult)) > 1
    fprintf('%s\n',strResult)
end
