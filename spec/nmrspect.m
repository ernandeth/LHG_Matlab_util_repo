function f = NMRspect( Mo, wo, T2 , w)
% function f = NMRspect(Mo, wo, T2 , w)
% This function returns the absorption of a spectrum made up
% of Lorentzian functions detrmined by the elements in the following arrays 
% magnitude (Mo)
% center frequency (wo)
% T2 
% and a range of frequencies (w)
% f = Mo*T2./(1 + ((w-wo).*T2).^2 );

Npeaks = size(Mo,2);
f=zeros(size(w));

for i=1:Npeaks
    f = f + Mo(i)*T2(i)./(1 + ((w-wo(i)).*T2(i)).^2 );
end

return
