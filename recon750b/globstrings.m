function results = globstrings(patterns, strs);
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: globstrings.m 604 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
%GLOBSTRINGS  String matching via wildcards.
% GLOBSTRINGS(PATTERNS,STRINGS) returns a cell array of the strings from 
% STRINGS that match some wildcard pattern in PATTERNS.
% STRINGS is a cell array of strings.
% PATTERNS is a string or cell array of strings.
%
% Two types of wildcards are supported:
% * matches zero or more characters.
% ? matches exactly one character.
%
% Examples:
%   globstrings('f?*r',{'fr','fur','four'})   % returns {'fur','four'}
%   globstrings({'a*','*c'},{'ace','bar','rac'}) % returns {'ace','rac'}
%
% See also glob, glob2regexp.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

if ischar(patterns)
  patterns = cellstr(patterns)';
end
patterns = patterns;
results = {};
for i = 1:length(patterns)
  s = regexp(strs,glob2regexp(patterns{i}),'match');
  s = cat(2,s{:});
  if isempty(s)
    warning([patterns{i} ' did not match anything']);
  else
    results = {results{:} s{:}};
  end
end
