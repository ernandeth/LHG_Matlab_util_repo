function fit_data = fit_spectra(data, T2, wo)

% function fit_data = fit_spectra(data, T2, wo)
%
% This function fits the amplitude and T2 of a 
% set of spectral peaks to the imput data
% the data must be in the form of a matrix where 
% each FID is a row.
verbose=0;
if (nargin==1)

    b=[0 21.44 85.7 193 343 536 772 1051]; %(s/mm2)

    %    water   mI    Cho   Cr    NAA
    T2 = [70     270   350   380   305]     /1000;
    wo = [4.7    3.59  3.21  3.04  2.02]; 
end

data = data';

% Account for scanner offset (not quite on 128 MHz)
shiftdrift = 0.16;

% Assuming Bandwith = 5000, Number of points/ FID = 2048
T2 = T2 * 5000/2048;

% Water resonance:
wo = [4.7 - wo + shiftdrift] * 128 * 2048/5000;
wo(1) = 0;

% Initialize some variables:
parms = [T2'  wo'];
fit_data = [];

for i=1:size(data,1)

    FID = data(i,:);
    
    % First remove the water peak by fitting its T2 decay
    % to the FID
	guess =lsqcurvefit('t2_decay', [16 0.001],[1:2048], abs(FID));
	m0 = guess(1);
	R2 = guess(2);
	water = t2_decay([m0 R2],  [1:2048]);
    
    dry_data = abs(FID)- abs(water);
    
    fdata=(abs(fft(dry_data)));
    fdata = fdata(1:1024);

    
    % Now do a fit of all the Lorentzians to each spectrum
    Mofit = lsqcurvefit('NMRfunc',[20 1 2 3 4 ], [1:1024], fdata, [],[],[],parms); 
    tmp = NMRspect(Mofit, wo, T2, [1:1024]);
    fit_data = [fit_data; tmp];
    if verbose==1
        clf
        plot(fdata)
        hold on
        plot(tmp, 'r');
        axis([50 170 0 1000])
        pause
    end
end

fit_data = fit_data';

return