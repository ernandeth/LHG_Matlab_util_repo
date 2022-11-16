function cooked= k_despiker(raw, skip, LEVEL, SHOWPIX)
%function [cooked,info] = k_despiker(pfile, skip, LEVEL, SHOWPIX)


fprintf('\nRunning despiker @ threshold : %f std devs', LEVEL);
%[raw,info] = read_raw_3d(pfile,0);
[nframes ndat ncoils] = size(raw);
cooked = raw;

for nc = 1:ncoils
    fprintf('\n\tDespiking coil %d', nc);
    
    odds  = squeeze(raw(skip+1:2:nframes, :,  nc));
    evens = squeeze(raw(skip+2:2:nframes, :,  nc));
    
    tmp = dspk(odds, LEVEL,SHOWPIX);
    cooked(skip+1:2:nframes, :, nc) = tmp;
    
    tmp = dspk(evens, LEVEL,SHOWPIX);
    cooked(skip+2:2:nframes, :, nc) = tmp;
    
end

return

function out = dspk(in, LEVEL, SHOWPIX)

out = in;
nfr = size(in,1);
m = mean(in, 1);
m = repmat(m,nfr,1);

s = std(in,[],1);
s = repmat(s,nfr,1);

% deviation from the mean
dev = abs(in-m);
% find those whos deviate more than LEVEL* standard deviations
inds = find(dev > LEVEL*(s));

for k=1: length(inds)
    [row col] = ind2sub( size(in), inds(k));
    
    if (row>1) & (row < nfr)
        out(row,col) = (in(row-1,col) + out(row+1,col) )/2;
        %fprintf('\n%d, %d- replacing row %d, col %d',k, inds(k), row, col);
    end
end
%
if(SHOWPIX)
    subplot(221)
    imagesc(abs(in))
    subplot(222)
    imagesc(abs(out))
    subplot(223)
    imagesc(abs(in-out))
    drawnow
end
%}
return