function pattern = glob2regexp(pattern)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: glob2regexp.m 603 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
%GLOB2REGEXP   Convert a glob pattern into a regexp pattern.
% GLOB2REGEXP(PATTERN) returns a regexp pattern which matches the same strings
% as the given glob pattern.
%
% Examples:
%   glob2regexp('my*file?.txt')
%   returns '^my.*file.\.txt$'

pattern = strrep(pattern,'.','\.');
pattern = strrep(pattern,'*','.*');
pattern = strrep(pattern,'?','.');
pattern = ['^' pattern '$'];
