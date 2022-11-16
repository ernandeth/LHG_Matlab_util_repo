function b_corr = phaseDetrend_220412(b,Ncenter,doSHOWPIX,fh)
%PHASEDETREND: function for correcting phase by evaluating drift at
%   center of kspace
%|
%|  Accepts "b-matrix" (raw data reshaped into {# data pts per echo} x
%|      {# of shots} x {# of slices}) and returns version of b corrected
%|      for phase drift
%|
%|  Parameters:
%|      b: "b-matrix" from spiralrecon
%|      Ncenter: number of points at center of kspace
%|      doSHOWPIX: option to display plots for debugging
%|      fh: figure handle to use
%|

    % Determine indicies of navigator points:
    ndat = size(b,1);
    navpts = (round(ndat/2 - Ncenter/2 + 1) : round(ndat/2 + Ncenter/2))';

    % Average together average phase at navigators and subtract it out:
    b_corr = abs(b).*exp(sqrt(-1)*(angle(b) - mean(angle(b(navpts,:,:)),1)));
    
    % Show how it affects phase:
    if doSHOWPIX
        figure(fh);
        plot(unwrap(angle(b(:)))), hold on
        plot(unwrap(angle(b_corr(:)))), hold off
        legend('Uncorrected','Corrected')
    end

return