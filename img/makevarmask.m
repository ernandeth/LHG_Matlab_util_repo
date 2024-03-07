function msk = makevarmask(raw , th)
% function msk = makevarmask(raw , th)
% 
% takes a 4D data set and 
% makes a mask that preserves the data with the top TH percentage of the
% variance in the data
% output: binary mask (3D)

dim=size(raw);
if numel(dim)==4
    varimage=std(raw,[],4);
elseif numel(dim)==3
    varimage=std(raw,[],3);
end

ordered = sort(varimage(:));
Nintgrl = cumsum(ordered)/sum(ordered(:)) * 100;
thval = find(Nintgrl>100-th);
%subplot(211), plot(ordered); title('Ordered Std. Deviations')
%subplot(212), plot(Nintgrl); title('Intrageted , Normalized Std. Deviations')

thval = ordered(thval(1))
msk = varimage;
msk(msk < thval) = 0;
msk(msk>=thval) = 1;

return