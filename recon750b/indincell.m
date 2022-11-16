function ind = indincell(cas,str)
% indincell.m - find index that a string is found in a cell array of strings
% INPUTS
% cas - cell array of strings
% str - string to find index of
% 
% OUTPUTS
% ind - integer or vector of integers of indices of matches
% 
% EXAMPLE
% cas = {'one';'three';'seven';'twelve'};
% indThree = indincell(cas,'three')

% Author - Krisanne Litinas
% $Id: indincell.m 573 2013-05-20 16:00:48Z klitinas $

iMatch = strfind(cas,str);
ind = findnonemptycells(iMatch);