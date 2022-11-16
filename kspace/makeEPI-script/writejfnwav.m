function writejfnwav(g, fname, hdr, scale, maxg)
% function writejfnwav(g, fname, hdr, scale, maxg)
%
% Use jfn_files.c to load .jfnwav files into your psd.
% See for example /server/home/jfnielse/epic/genericpsd/pc.e.
%
% INPUTS:
%   g      - gradient waveform (G/cm)
%   fname  - output .jfnwav file name
%   hdr    -  hdr.res_g
%             hdr.res_kpre
%             hdr.res_k
%             hdr.npix
%             hdr.nechoes
%   maxg    - max gradient (G/cm - usually 4)
%
% jfn, 2006
% $Id: writejfnwav.m,v 1.1 2010/03/15 02:18:35 jfnielse Exp $


% scale waveform to +- 32766
g = g/maxg;

g = 2*round(scale*g*32766/2); 

% set EOS bit of last point
if g(length(g)) == 0
	g(length(g)) = 1;
else
	g(length(g)) = g(length(g)) - sign(g(length(g)));
end

fid = fopen(fname, 'w');
fwrite(fid, hdr.res_g,    'int16', 'ieee-be');
fwrite(fid, hdr.res_kpre, 'int16', 'ieee-be');
fwrite(fid, hdr.res_k,    'int16', 'ieee-be');
fwrite(fid, hdr.npix,     'int16', 'ieee-be');
fwrite(fid, hdr.nechoes,  'int16', 'ieee-be');
fwrite(fid, hdr.nechosep, 'int16', 'ieee-be');
fwrite(fid, g,            'int16', 'ieee-be');
fclose(fid);

return;

