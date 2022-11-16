function casFile = getverifyfiles(strPattern,strDir,strType,rng)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: getverifyfiles.m 601 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strPattern = '.*_0{1,2}[1,2]\.img$';
% strDir = 'C:\temp\fmri_testing\c1317plas_control\study_20050909\series_10_bas_MoCoSeries';
% strType = 'cas';
% rng=1:4; % from one thru four files expected
% casFiles = getverifyfiles(strPattern,strDir,strType,rng)

casFile = getpatternfiles(strPattern,strDir,strType);
if ~verifyfilecount(casFile,strType,rng)
    error('unexpected file count')
end
