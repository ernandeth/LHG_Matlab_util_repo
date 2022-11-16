function [D , inds] = mk_overlay(threshold, D, spm_data, minBin)

% BUT if fthe threshold is negative, we threshold in the NEGATIVE DIRECTION
    if (threshold < 0)
        spm_data = -spm_data;
        fprintf('\n WARNING: SPM NEGATIVE THRESHOLD -> WILL SELECT MOST NEGATIVE VALUES\n');
    end
%     inds = find(spm_data >= abs(threshold));
inds = find(spm_data > abs(threshold));
    
    % window the map
%     Bmax = max(abs(spm_data(inds)));
    Bmax = max(spm_data(inds));

    Bmin = abs(threshold);
    range = Bmax-Bmin;
    if isempty(range)
        range=1;
	end
	
    % scale the activation maps to use the whole colormap
    spm_data = (spm_data - Bmin )* 256 /range ;
    spm_data(find(spm_data < abs(threshold)))=0;


    %spm_data2(inds) = 0;

    if ~isempty(inds)
        D(inds) =  minBin + spm_data(inds);
    end