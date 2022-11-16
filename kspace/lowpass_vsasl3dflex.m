function raw_filt = lowpass_vsasl3dflex(raw,echoLowpass)

% Get dimensions of raw
nframes = size(raw,1);
ndat = size(raw,2);
nleaves = size(raw,3);
nslices = size(raw,4);
ncoils = size(raw,5);

% Create polynomial coefficients for filter (code format copied from
% filterDesigner generated code)
[N, Fo, Ao, W] = firpmord(echoLowpass*[0.9 1.05], [1 0], [5e-3, 1e-4]);
b = firpm(N, Fo, Ao, W, 20);
a = 1;

% The following code is a modified version of filtfilt(), but acts on the second
% dimension of data to be compatible with raw's format

        len = ndat;
        b = b(:).';
        a = a(:).';
        nb = length(b);
        na = length(a);
        nfilt = max(nb,na);

        nfact = 3*(nfilt-1);  % length of edge transients
        if (len<=nfact)    % input data too short!
            error('\nRaw data must have length more than 3 times filter order. Try a higher echoLowpass');
        end

    % set up filter's initial conditions to remove dc offset problems at the 
    % beginning and end of the sequence
        if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
        if na < nfilt, a(nfilt)=0; end

    % use sparse matrix to solve system of linear equations for initial conditions
    % zi are the steady-state states of the filter b(z)/a(z) in the state-space 
    % implementation of the 'filter' command.
        rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
        cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
        data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
        sp = sparse(rows,cols,data);
        zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

    % Extrapolate beginning and end of data sequence using a "reflection
    % method".  Slopes of original and extrapolated sequences match at
    % the end points.
    % This reduces end effects.
        raw_filt = cat(2, ...
            2*raw(:,1,:,:,:)-raw(:,(nfact+1):-1:2,:,:,:), ...
            raw, ...
            2*raw(:,len,:,:,:)-raw(:,(len-1):-1:len-nfact,:,:,:));

    % filter, reverse data, filter again, and reverse data again
        raw_filt = filter(b,a,raw_filt,zi*raw_filt(1),2);
        raw_filt = raw_filt(:,length(raw_filt):-1:1,:,:,:);
        raw_filt = filter(b,a,raw_filt,zi*raw_filt(1),2);
        raw_filt = raw_filt(:,length(raw_filt):-1:1,:,:,:);

    % remove extrapolated pieces of y
        raw_filt(:,[1:nfact len+nfact+(1:nfact)],:,:,:) = [];
    
    
end

