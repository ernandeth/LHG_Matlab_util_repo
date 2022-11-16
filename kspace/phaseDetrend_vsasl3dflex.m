function [S_corr,fits] = phaseDetrend_vsasl3dflex(S,navpts,order,isSOS)
% function [S_corr,fits] = phaseDetrend_vsasl3dflex(S,navpts,order,isSOS)
%
% S :         the input data block
% navpts ;    which points from each echo are being used for correction,
%             usually  center of k space
% order:      order of the polynomial to be regrressed out
% isSOS:      indicate whether it's a stack of spirals
%
%   Data dimensions
% 	nframes = size(S,1);
% 	ndat = size(S,2);
% 	nleaves = size(S,3);
% 	nslices = size(S,4);
% 	ncoils = size(S,5);

    % Extract dimensions
	nframes = size(S,1);
	ndat = size(S,2);
	nleaves = size(S,3);
	nslices = size(S,4);
	ncoils = size(S,5);
    
    % Initialize S_corr and fits
    S_corr = zeros(size(S));
    fits = S_corr;
    
    % Create design matrix from navpts and endpoints, centered around
    % middle of echo
    if ~isSOS % if SERIOS, include endpoints as navigators
%         navpts = [1;navpts(:);ndat];
    end
    x_nav = navpts  - round(ndat/2);
    x_all = (1:ndat) - round(ndat/2);
    A = x_nav(:).^(order:-1:0);
    
    % Loop through each echo and fit a polynomial to the phase drift at
    % navigators
    for framen = 1:nframes
        for leafn = 1:nleaves
            for slicen = 1:nslices
                for coiln = 1:ncoils
                    echo = S(framen,:,leafn,slicen,coiln);
                    y = unwrap(angle(echo(navpts)));
                    betas = pinv(A)*y(:);
                    fits(framen,:,leafn,slicen,coiln) = ...
                        x_all(:).^(order:-1:0) * betas;
                end
            end
        end
    end
    
    % Correct signal; if SOS, add back average phase at center of echo for each
    % slice/leaf
    S_corr = S.*exp(-1i*fits + 1i*isSOS*mean(angle(S(:,navpts,:,:,:)),[1,2,5]));
    
end
