function t = gent(v,fs,offset)
% General-Purpose Toolbox - a package of utilities used primarily for
% in-house analysis workflows.  This set of tools was developed as an
% in-house package by the Stroke Rehab Team of the Functional Electrical
% Stimulation (FES) Center in Cleveland, Ohio. This software is available
% to the scientific community as copyright freeware under the terms of the
% GNU General Public License (GPL).
%
% $Id: gent.m 624 2013-05-20 16:38:02Z klitinas $
%
% Authors: Sahil Grover, Kenneth Hrovat, Krisanne Litinas, Roger Cheng (c)
% Stroke Rehab Team of the Functional Electrical Stimulation (FES) Center
% in Cleveland, Ohio.  All rights reserved.
% GENT generate time vector, t, (in seconds) to match length of input vector & sample rate (in samples per second) 
%
% USAGE:
% t = gent(v,fs[,offset]);
%
% INPUTS:
% v - vector of values to create time vector for based on its length
% fs - scalar sample rate in samples/second
% offset - optional scalar time offset (in seconds) such that t starts at this value (and not at default of zero)
%
% OUTPUTS:
% t - vector of time values (in seconds) based on fs and length of v
%
% EXAMPLE:
% t1 = gent(1:5,10)
% t2 = gent(1:5,10,95)

% Author: Ken Hrovat
% $Id: gent.m 624 2013-05-20 16:38:02Z klitinas $

% Generate time vector starting at zero with step size (delta t) from fs & duration established by length of v
if nargin > 3 || nargin < 2
    error('daly:common:wrongNumberOfArguments','incorrect number of inputs supplied to %s',mfilename)
else
    t = 0:1/fs:length(v)/fs-1/fs; % in seconds
    % Offset time vector (if needed)
    if nargin == 3
        t = t + offset; % in seconds
    end
    % As courtesy to calling routine, reshape t to match input vector v
    t = reshape(t,size(v)); % in seconds
end
