function [bgst1 bgst2] = optimal_bgs_time(tdelay, t1)
% calculate BGS inversion times using the equation from Guenther 2005 -
% MRM54(2):491-8
if nargin==0
    tdelay = 1.5
    t1 = 0.7
end

bgst1 = tdelay + 2*t1*log( (1/2 + 1/4) + (1/2 - 1/4)*exp(-t1*tdelay/2))
bgst2 = tdelay + 2*t1*log( (1/2 - 1/4) + (1/2 + 1/4)*exp(-t1*tdelay/2))

fprintf('\nBGS at %f and %f before readout \nwill suppress tissues with T1=%f and %f \nat TI=%f',...
    bgst1, bgst2, t1, 2*t1, tdelay)
fprintf('\nBS1_time = %f \nBS2time = %f', tdelay-bgst1, bgst1-bgst2);
return