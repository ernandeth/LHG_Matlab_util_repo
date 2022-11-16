function [g] = makecrusher(cycles, opslthick)
% function [g] = makecrusher(cycles, opslthick)
%
% Generate crusher.pcwav
%
%  opslthick  - mm
%
% Use with ~jfnielse/genericpsd/pc.e
%
% $Id: makecrusher.m,v 1.1 2010/03/15 02:18:35 jfnielse Exp $

opslthick = opslthick/10;              % cm
gamma = 4257.5;                        % Hz/Gauss
obl = sqrt(4);
mxg = 3.9/obl;                         % Gauss/cm
mxs = 14.5/obl;                        % Gauss/cm/msec
t = cycles/(opslthick*gamma*mxg)*1e3;  % msec
n = round(t/4e-3);                     % assumes 4 us sampling

s = mxs * 4e-3;                        % max change in g per sample (G/cm)
g = [s:s:mxg mxg*ones(1,n) (mxg-s):-s:0];
g = [g zeros(1,mod(length(g),2))];

% write wav files
hdr.res_g = length(g);
hdr.res_kpre = 0;
hdr.res_k = 0;
hdr.npix = 0;  
hdr.nechoes = 0;  
hdr.nechosep = 0;

prefix = [num2str(cycles) 'cycles_' num2str(opslthick) 'slthick_crusher'];

% writejfnwav(g, [prefix '.jfnwav'], hdr, 1.0, 4.0)
writejfnwav(g, ['crusher.jfnwav'], hdr, 1.0, 4.0)

plot(g);

% EOF
