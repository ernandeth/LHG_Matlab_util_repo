function [N,S]=nec_suf(x,y)
% function [N,S]=nec_suf(x,y)
%
% Uses formulas from Judea Pearl's causality book.

t = [sum(x&y) , sum(~x&y);
    sum(x&~y) , sum(~x&~y)];

sumrow = sum(t);

p_y_given_x = t(1,1)/sumrow(1);

p_y_given_notx = t(1,2)/sumrow(2);

p_noty_given_notx = t(2,2)/sumrow(2);

PNS = p_y_given_x - p_y_given_notx;

bound1 = max(0 , p_y_given_x - p_y_given_notx);

bound2 = min(p_y_given_x , p_noty_given_notx);

if PNS<bound1
    PNS=0;
    N=0;
    S=0;
    return
    
elseif PNS>bound2
    PNS=1;
    N=1;
    S=1;
    return
end

N = PNS / p_y_given_x;
S = PNS / (1 - p_y_given_notx);

if isnan(N) | isinf(N)
    N=0;
end

if isnan(S) | isinf(S)
    S=0;
end
