function casFiles = recursivelyfindfiles(strTop,strPattern)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: recursivelyfindfiles.m 614 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strTop = 'C:\data\fmri\fromopticalmedia\temp4preone\s1372plas\13720000';
% strPattern = '.*\.dcm';
% casFiles = recursivelyfindfiles(strTop,strPattern);

%% Get parameters to work with spm_select
sFilt.code = 0;
sFilt.frames = [];
sFilt.ext = {'.*'};
sFilt.filt = {strPattern};

%% Call recurse functionality added to spm_select
casFiles = spm_select('recurse',strTop,sFilt);
