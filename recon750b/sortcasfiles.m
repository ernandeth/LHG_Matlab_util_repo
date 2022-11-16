function [casSort,iSort] = sortcasfiles(casFiles,strPrefix)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: sortcasfiles.m 618 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% SORTCASFILES.M - takes care of files with unsorted numbering or
% different length of numbers [i.e., 99 and 100 instead of 099 and 100]
% 
% INPUTS
% casFiles - list of files to sort
% 
% OUTPUTS
% casSort - list of sorted filenames [same size as casFiles]
% 
% EXAMPLE
% strDir = 'S:\data\upper\vicon\dalyUE\upperStroke\s1370plas\20091221_s1370plas';
% strPat = 'Supination_pronation.*c3d$';
% casFiles = getpatternfiles(strPat,strDir,'cas');
% strPrefix = 'supination_pronation';
% casSort = sortcasfiles(casFiles,strPrefix);

% Author - Krisanne Litinas
% $Id: sortcasfiles.m 618 2013-05-20 16:16:56Z klitinas $

% convert strPrefix to cell, then repeat to match size of casFiles
casPrefix = repmat({strPrefix},length(casFiles),1);

% call locgetnum to do regular expression matching
cas = cellfun(@locgetnum,casFiles,casPrefix,'uni',0);

% str2double the cas
iFile = cellfun(@str2double,cas);

[foo, iSort] = sort(iFile);
casSort = casFiles(iSort);

%-----------------------------------
function n = locgetnum(str,strPrefix)
strBase = basename(str);
strPattern = ['^' strPrefix '\D*(?<num>\d*)$'];
x = regexpi(strBase,strPattern,'names');
n = x.num;
