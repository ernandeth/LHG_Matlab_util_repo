function fileSize = joliet(dirName,nLength)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: joliet.m 606 2013-05-20 16:16:56Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
%JOLIET Check a directory prior to burning a CD
%   FILESIZE = JOLIET(DIRNAME,NAMELIMIT)
%   FILESIZE = JOLIET(NAMELIMIT)
%   FILESIZE = JOLIET(DIRNAME)
%
%   The "Joliet" filesystem used by CD's limits file length to 64 characters.
%   This m-file will traverse a directory and all sub-directories and warn if
%   any files have names longer than 64 characters.  If DIRNAME is omitted, 
%   the function will use the present working directory (PWD).  The optional
%   NAMELIMIT argument sets a length other than 64 for warning.  The function
%   returns the total size of all files in megabytes as FILESIZE.  Most CD's
%   can contain 600-700 megabytes of data.
%
%   In the following example, we see that there is a track with a 73 character
%   long name in the user's audio directory.  This should be shortened before
%   buring a CD.  The directory contains enough data to fill a CD halfway.
%
%   >> joliet c:\audio
%         > c:\audio
%         > c:\audio\Nine Inch Nails
%         > c:\audio\Nine Inch Nails\Halo 03 - Head Like a Hole
%         > c:\audio\Nine Inch Nails\Halo 04 - Sin
%         > c:\audio\Nine Inch Nails\Halo 06 - Fixed
%         > c:\audio\Nine Inch Nails\Halo 07 - March of the Pigs
%   (73) Nine Inch Nails - March of the Pigs - 03 - All the Pigs, All Lined Up.mp3
%         > c:\audio\Nine Inch Nails\Halo 09 - Closer to God
%         > c:\audio\Nine Inch Nails\Halo 15 - We're In This Together
%    ans = 
%       313.2
%
%   @AUTHOR  Karl Critz, The MathWorks
%   @VERSION 1.0

defaultDir = pwd;
if nargin<2
    nLength = 64;
    if nargin<1
        dirName = defaultDir;
    elseif isnumeric(dirName)
        nLength = dirName;
        dirName = defaultDir;
    elseif ~isnan(str2double(dirName))
        nLength = str2double(dirName);
        dirName = defaultDir;
    end
end

%disp(['        > ' dirName]);

allFiles = dir(dirName);

fileIdx = find(~[allFiles.isdir]);
for i=fileIdx
    if length(allFiles(i).name)>nLength
        disp(['        > ' dirName]);
        disp(sprintf('(%i) %s',length(allFiles(i).name),allFiles(i).name));
    end
end
fileSize = sum([allFiles(fileIdx).bytes])/1e6;

dirIdx = find([allFiles.isdir]);
for i=dirIdx
    if ~any(strcmp({'.','..'},allFiles(i).name));
        if length(allFiles(i).name)>nLength
            disp(sprintf('(%i) %s',length(allFiles(i).name),allFiles(i).name));
        end
        fileSize = fileSize + joliet(fullfile(dirName,allFiles(i).name),nLength);   
    end
end
%format short g
