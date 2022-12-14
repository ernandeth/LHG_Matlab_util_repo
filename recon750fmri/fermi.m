function f= fermi(sz,radius,width)
%  usage ...  fermi(sz,radius,width)
%  sz - matrix size, radius - radius of window, width - transition width
%  negative width will give "standard width" = 20*radius/128

%
% [x,y]= meshdom(-64:1:63,  63:-1:-64);
% f= 1 ./ (1 + exp((sqrt(x.^2 + y.^2) - radius) ./ (10*steepness/256)));

% $Id: fermi.m 1266 2014-03-20 18:12:40Z klitinas $
% $HeadURL: svn+ssh://klitinas@woo.engin.umich.edu/work/repos/matlab/recon/branches/phase_zero/fermi.m $

if width < 0,
   width = 20*radius/128;
end
i = sqrt(-1);
cent = sz/2 + 1;
x= (1:sz);
y= (1:sz)';
X= ones(size(y))*x;
Y= y*ones(size(x));
clear x y;
R =  abs( X-cent + (Y-cent).*i );
clear X Y;
f = 1 ./ (1 + exp( (R - radius) ./ width ));
clear R;
