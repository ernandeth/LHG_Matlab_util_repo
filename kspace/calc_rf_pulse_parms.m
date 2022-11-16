function [paramsfloat rfstruct] = sub_myrfstat(b1, nom_fa, system_maxRF);
% Calculate RF parameters needed for RFPULSE struct in .e file.
% Needed for B1 scaling, SAR calculations, and enforcing duty cycle limits.
% See also mat2signa_krishna.m
%
% b1         real 1D vector containing B1 amplitude, size Nx1 [Gauss]
% nom_fa     nominal flip angle (degrees)
%
% RF waveform is scaled relative to system.maxRf.
% This may be 0.25G/0.125G for quadrature/body RF coils (according to John Pauly RF class notes), but haven't verified...

g = 1;  % legacy dummy value, ignore
nom_bw = 2000;

dt = 4e-6;                        % use 4 us RF sample width
%gamma = 4.2575e3;                  % Hz/Gauss
tbwdummy = 2;

hardpulse = max(abs(b1)) * ones(length(b1),1);    % hard pulse of equal duration and amplitude

pw            = length(b1)*dt*1e3;                                       % ms
nom_pw        = length(b1)*dt*1e6;                                        % us

if max(abs(b1)) == 0   % non-zero RF pulse
	error('RF waveform cannot be zero');
end
abswidth      = sum(abs(b1)) / sum(abs(hardpulse));
effwidth      = sum(b1.^2)   / sum(hardpulse.^2);
% or equivalently:  effwidth = sum(b1.^2)/(max(abs(b1))^2)/length(b1)
area          = abs(sum(b1)) / abs(sum(hardpulse)); 
dtycyc        = length(find(abs(b1)>0.2236*max(abs(b1)))) / length(b1);
maxpw         = dtycyc;
num           = 1;
max_b1        = system_maxRF; % Gauss. Full instruction amplitude (32766) should produce max_b1 RF amplitude,
					% as long as other RF .mod files (if any) use the same max_b1.
max_int_b1_sq = max( cumsum(abs(b1).^2)*dt*1e3 );   	% Gauss^2 - ms
max_rms_b1    = sqrt(mean(abs(b1).^2));              	% Gauss
nom_fa        = nom_fa;                              	% degrees
%nom_bw        = tbwdummy / (dt * length(b1));       	% Hz
% max_int_b1    = abs(sum(b1))*dt*1000

% calculate equivalent number of standard pulses
stdpw = 1;                                                % duration of standard pulse (ms)
stdpulse = 0.117 * ones(round(stdpw/(dt*1e3)),1);
numstdpulses = num * effwidth * (pw/stdpw) * (max(abs(b1))/0.117)^2;

%pulse50 = 0.117/180*50 * ones(round(stdpw/(dt*1e3)),1);
%num50 = sum(abs(b1).^2)/sum(pulse50.^2);

paramsfloat = [pw     abswidth  effwidth area          dtycyc      ...
              maxpw  num       max_b1   max_int_b1_sq max_rms_b1  ...
              90 nom_pw    nom_bw   g numstdpulses  nom_fa            ];  % hardcode 'opflip' to 90

rfstruct.pw = pw;
rfstruct.abswidth = abswidth;
rfstruct.effwidth = effwidth;
rfstruct.area = area;
rfstruct.effwidth = effwidth;
rfstruct.dtycyc = dtycyc;
rfstruct.maxpw = maxpw;
rfstruct.num = num;
rfstruct.max_b1 = max_b1;
rfstruct.max_int_b1_sq = max_int_b1_sq;
rfstruct.max_rms_b1 = max_rms_b1;
rfstruct.nom_90 = 90;
rfstruct.nom_pw = nom_pw;
rfstruct.nom_bw = nom_bw;
rfstruct.g = g;
rfstruct.numstdpulses = numstdpulses;
rfstruct.nom_fa = nom_fa;

rfstruct

return