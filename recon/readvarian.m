% %   Function to upload an fid recorded by vnmrj
%        [Fid, procpar] = getFid(directory)
%         Inputs:  directory:  string location of the fid and it's associated parameter file
%                   
%         Outputs: Fid:  Free induction decay.  a nxm matrix of complex
%                            points that is the free induction decay recorded by vnmrj
%                  procpar:  Structure containing Vnmrj acquisition parameters
function [Fid, procpar] = getFid(directory)

procpar = tryreadprocpar(directory);
[Fid, ~, ~] = readfid([directory '/fid'], -1);

end