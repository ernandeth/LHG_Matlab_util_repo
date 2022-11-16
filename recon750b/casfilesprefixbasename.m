function casOut = casfilesprefixbasename(casFiles,strPrefix)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: casfilesprefixbasename.m 579 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% casFiles = {[];'c:\temp\swuaOne.img';'c:\temp\swuaTwo.img'};
% casOut = casfilesprefixbasename(casFiles,'d');

casOut = casFiles;
for i = 1:length(casFiles)
    strFile = casFiles{i};
    if isempty(strFile), continue; end
    [strPath,strName,strExt] = fileparts(strFile);
    strBase = [strPrefix strName strExt];
    casOut{i} = fullfile(strPath,strBase);
end
