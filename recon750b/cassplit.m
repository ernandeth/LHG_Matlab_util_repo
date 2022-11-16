function casNew = cassplit(cas, strSplitter)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: cassplit.m 571 2013-05-20 15:59:43Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% cassplit - splits each element of cell array of strings using a splitter
% 
% INPUTS
% cas - nx1 cell array of strings to be split
% strSplitter - string, where you want the rows to be split
% 
% OUTPUTS
% casNew - nxm cell array of strings, where m is the result of split
% 
% NOTE
% I'm in a hurry, so this only works if there is the same number of
% splitters in each row.
% 
% EXAMPLE
% cas = {'one\two\three'; 'four\five\six'};
% casNew = cassplit(cas,'\')

% Author - Krisanne Litinas
% $Id: cassplit.m 571 2013-05-20 15:59:43Z klitinas $

% Split and unnest the cellfun-determined result
casSplitter = repmat({strSplitter},size(cas));
casNest = cellfun(@strsplit,casSplitter,cas,'uni',false);
casUnnested = unnestcell(casNest);

% Determine size of output cas
numOutputRows = nRows(cas);
numOutputCols = numel(casUnnested)/numOutputRows;

% Reshape and transpose [kind of clumsy]
casReshaped = reshape(casUnnested,numOutputCols,numOutputRows);
casNew = casReshaped';
