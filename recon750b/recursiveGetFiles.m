function [files, dirs] = recursiveGetFiles(rootPath, varargin)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: recursiveGetFiles.m 613 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% RECURSIVEGETFILES -- function to recursively return filepaths and directories from a root path
% 
% Inputs:
%     rootPath -- The top level directory to search
%     extension(OPTIONAL) -- file extension filter
% 
% Outputs:
%     files -- cell array of strings of filepaths in top level directory and all sub directories
%     dirs -- all sub directories in top level directory
%     
% Example:
%     datadir = 'S:\data\upper\bci\therapy\c9999rand';
%     [files, dirs] = recursiveGetFiles(datadir, '.prm');

% Author: Sahil Grover

% optional file extension filter
ext = '*';
if nargin>1
    ext = [varargin{1}];
    if isempty(strfind(ext,'*'))
        ext = ['*' ext];
    end
end

[filenames, details] = dirdeal(rootPath);

% get indices of directories
[dirInd{1:length(details)}] = deal(details.isdir);
dirInd = cell2mat(dirInd);
filenames = filenames(dirInd);
dirs =  cellfun(@(x)[rootPath filesep x],filenames,'uni',0);

[rawFiles,fDetails] = dirdeal([rootPath filesep ext]);
[filesInd{1:length(fDetails)}] = deal(fDetails.isdir);
filesInd = ~cell2mat(filesInd);
files = rawFiles(filesInd);
files = cellfun(@(x)[rootPath filesep x],files,'uni',0);

for i=1:length(dirs)
    [f, d] = recursiveGetFiles(dirs{i},ext);
    files = [files; f];
    dirs = [dirs; d];
end
