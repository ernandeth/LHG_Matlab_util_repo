function [g, hdr] = readjfnwav(fname)
% function [g, hdr] = readjfnwav(fname)
%
% Use jfn_files.c to load .jfnwav files into your psd.
%
% jfn, 2006
% $Id: readjfnwav.m,v 1.1 2010/03/15 02:18:35 jfnielse Exp $

fid = fopen(fname, 'r');
hdr.res_g    = fread(fid, 1,  'int16', 0, 'ieee-be');
hdr.res_kpre = fread(fid, 1,  'int16', 0, 'ieee-be');
hdr.res_k    = fread(fid, 1,  'int16', 0, 'ieee-be');
hdr.npix     = fread(fid, 1,  'int16', 0, 'ieee-be');
hdr.nechoes  = fread(fid, 1,  'int16', 0, 'ieee-be');
hdr.nechosep = fread(fid, 1,  'int16', 0, 'ieee-be');
g = fread(fid, hdr.res_g,  'int16', 0, 'ieee-be');
fclose(fid);

return;

