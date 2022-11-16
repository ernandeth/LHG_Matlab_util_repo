function mymap = rt_make_colormap

mygray = [0:256]' * [1 1 1];

myhot = [128.5:0.5:255]' * [1 0 0] ;
tmp =   [128.5:0.5:255]' * [0 1 0] ;
tmp(:,1) = 256;
myhot = [myhot; tmp];
tmp =   [128.5:0.5:255]' * [0 0 1];
tmp(:,1:2) = 256;
myhot =  [myhot;  tmp;];
myhot = myhot(1:3:end, :);

myblue = myhot;
myblue(:,1) = myhot(:,3);
myblue(:,3) = myhot(:,1);

mymap = [mygray; myhot; myblue]/256;

Ncolors = size(myhot,1);
% if isstr(root)
%       set(gcf, 'Name', root)
% end

return
