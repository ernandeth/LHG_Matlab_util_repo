function f = NMRfunc(guess, w, parms)
% function f = NMRfunc(guess, w, parms)
% simulates an FID based on the given parameters for a set 
% of nuclei

Mo=guess;
T2 = parms(:,1);
wo = parms(:,2);

f = NMRspect(Mo,wo,T2,w);

return
