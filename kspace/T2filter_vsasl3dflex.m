function [S_corr,T2,T2curve] = T2filter_vsasl3dflex(S,T2weight,navpts,TE,dt)
% function [S_corr,T2,T2curve] = T2filter_vsasl3dflex(S,T2weight,navpts,TE,dt)
%
% Fits the echo train to a T2 decay curve and then multiplies the echo
% train by a compensation function to correct for the decay , like this:
% echo_out(n) = echo(n)*exp(R2*T2weights* n)
%
% It uses only the center of k-space (navpoints) to fit the T2 decay curve
%
% S:            input data block
% T2weight:     weight of the compensation f(n)
% navpts:       which points we use for fitting the T2 decay curce
% T2:           echo time
% dt:           time resolution for T2 decay curve
%

	% Initialize S_corr
    S_corr = zeros(size(S));
    % Extract dimensions
	nframes = size(S,1);
	ndat = size(S,2);
	nleaves = size(S,3);
	nslices = size(S,4);
	ncoils = size(S,5);
	ndat_all = round((TE*1e-3)/dt);
	
	% Initialize array of R2 values
	R2s = zeros(nframes,nleaves,ncoils);

	% Compute navigator signal for each interleaf/slice/frame/coil by
	% averaging the FID over navigator (center) points
	S_nav = mean(S(:,navpts(:),:,:,:),2);
    
    % Create design matrix using navigator locations
	x_nav = round(ndat/2) + ndat_all * (0:nslices-1)';
    A = x_nav.^[1 0];
    
    % Loop through all echo trains
    for framen = 1:nframes
        for leafn = 1:nleaves
            for coiln = 1:ncoils
                
                % Determine signal at navigator points for current echo
                y_nav = abs(reshape(S_nav(framen,:,leafn,:,coiln),[],1));
                
                % Fit exponential decay curve using least squares
                b = (A'*A)^(-1) * A' * log(y_nav);
                
                % Save R2 value to buffer
                R2s(framen,leafn,coiln) = b(1);
                
            end
        end
    end
    
    % Average together all R2 values to determine single value
    R2 = mean(R2s(:));
    T2 = dt*1e3/abs(R2);
    
    % Create a correction array with weighting
    T2curve = reshape(exp(x_nav*R2),[1,1,1,nslices,1]);
    
    % Loop through all echo trains
    for framen = 1:nframes
        for leafn = 1:nleaves
            for coiln = 1:ncoils
                
                % Correct signal by dividing by the T2 curve
                S_corr(framen,:,leafn,:,coiln) = ...
                    (1-T2weight) * S(framen,:,leafn,:,coiln) + ...
                    T2weight * S(framen,:,leafn,:,coiln) ./ ...
                    T2curve(1,1,1,:,1);
                
                % Normalize energy of corrected signal
                S_corr(framen,:,leafn,:,coiln) = ...
                    S_corr(framen,:,leafn,:,coiln) * ...
                    norm(reshape(S(framen,:,leafn,:,coiln),[],1),1) / ...
                    norm(reshape(S_corr(framen,:,leafn,:,coiln),[],1),1);
                
            end
        end
    end

    % Format T2curve like raw for easy plotting later
    T2curve = exp(b(2))*repmat(T2curve,[nframes,ndat,nleaves,1,ncoils]);

end
