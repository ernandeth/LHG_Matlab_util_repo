function [s,ernst] = ernst_eq(tr,t1,flip)
% function [s,ernst] = ernst_eq(tr,t1,flip)

flip = flip/180*pi;
E1 = exp(-tr./t1);
s = (1-E1).*sin(flip) ./ (1-cos(flip).*E1);

ernst = acos(E1)/pi*180;
return