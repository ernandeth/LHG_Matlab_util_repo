function zth = getZthreshold(Zname, rate)      
% function zth = getZthreshold(Zname, rate)      

[Zimg  h]   = read_img(Zname);
      % read t image, with df degrees of freedom
      Zimg(~isfinite(Zimg(:))) = [];

      P        = 1-cdf('norm',Zimg(:),0,1);

      [pID pN] = FDR(P,  rate);

      zID      = icdf('norm',1-pID,0,1);     % z threshold, indep or pos. correl.
      tN       = icdf('norm',1-pN,0,1);      % z threshold, no correl. assumptions
      
      return
      
%%
      function [pID,pN] = FDR(p,q)
% FORMAT [pID,pN] = FDR(p,q)
% 
% p   - vector of p-values
% q   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
%______________________________________________________________________________
% $Id: FDR.m,v 1.1 2009/10/20 09:04:30 nichols Exp $


p = p(isfinite(p));  % Toss NaN's
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
pN = p(max(find(p<=I/V*q/cVN)));
return