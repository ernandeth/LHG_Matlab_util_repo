function casFiles = findregex(strDir,strRegex)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: findregex.m 586 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% EXAMPLE
% strDir = 'c:\data\fmri\adat';
% strSession = 'preone';
% strTask = 'wrist';
% strRegex = ['.*[cs][0-9][0-9][0-9][0-9][a-z][a-z][a-z][a-z].' strSession '.study_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].series_[0-9][0-9]_' strTask '_.*_series_structure_' strTask '\.mat$'];
% casFiles = findregex(strDir,strRegex)

%% Get cygpath for use with find
strCmdCygpath = ['C:\cygwin\bin\cygpath.exe "' strDir '"'];
[foo,strDirCyg] = dos(strCmdCygpath);
strDirCygpath = deblank(strDirCyg);

%% Do find on mat files (returning cygpaths)
% find /cygdrive/c/data/fmri/adat -type f -regex '.*[cs][0-9][0-9][0-9][0-9][a-z][a-z][a-z][a-z].pre.study_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].series_[0-9][0-9]_shoulder_.*_series_structure_shoulder\.mat$' -exec cygpath -w {} \;
strCmd = ['C:\cygwin\bin\find.exe ' strDirCygpath ' -type f -iregex ' strRegex];
fprintf('\nUsing following cygwin command to find files:\n%s...',strCmd) 
[foo,strBig] = dos(strCmd);

%% Split str output from dos call into cas
casFiles = strsplit(char(10),strBig); % char(10) is CR/LF

%% Remove empty cells
casFilesCygpath = casFiles(findnonemptycells(casFiles))';

%% Convert to cygpath
casFiles = cellfun(@cygpath,casFilesCygpath,'uni',false);
fprintf('\ndone.\n')
