function [result_struct, bopen] = findreplace(f, srchtxt, rplctxt)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: findreplace.m 587 2013-05-20 16:16:55Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% Function FINDREPLACE
%           modifies a text file, replacing srchtxt with rplctxt
%
% Inputs:   f - (string) filename or directory name to search
%           srchtxt (string) string text to find
%           rplctxt (string) string replacement string
%
% Outputs:  result_struct (structure) representing search results, includes
%                                     all original lines of file, and
%                                     line_numbers contains the array of
%                                     line numbers where occurrences were
%                                     found.
%           bopen (integer) representing success or failure of fopen() call
%
% Auth: Matthias Beebe
%
% Last Modified:  3/07
%
if(nargin == 2)
    rplctxt = '';
    no_replace = 1;
else
    if(rplctxt == -1)
        no_replace = 1;
    else
        no_replace = 0;
    end
end

result_struct = struct('lines', [], 'line_nums', []);

bopen = 0;

try
    if(no_replace)
        fid = fopen(f, 'r');
    else
        fid = fopen(f, 'r+');
    end
    if(fid == -1)
        bopen = -1;
        return;
    end
catch
    return;
end
    
lines = {};

while 1
    tline = fgets(fid);
    if ~ischar(tline),  break, end
    lines = [lines; tline];
end

search_results = regexp(lines, srchtxt);

for i = 1:length(search_results)
    if(isempty(search_results{i}) == 0) %search_results is not empty
        result_struct(1).line_nums = [result_struct.line_nums; i];
    end
end

% fix all strings where replacement is needed
if(size(result_struct) > 0)
    if(no_replace ~= 1)
        for i = 1:length(result_struct.line_nums)
            str = regexprep(lines{result_struct.line_nums(i)}, srchtxt, rplctxt);
            lines{result_struct.line_nums(i)} = str;
        end

        % set file pointer to beginning of file
        frewind(fid);

        for i = 1:length(lines);
            fprintf(fid, '%s', lines{i});
        end

    end
end

result_struct(1).lines = lines;

fclose(fid);

end

