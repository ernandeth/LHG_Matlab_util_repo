function set_func_colors

%
% (c) 2005 Luis Hernandez-Garcia
% University of Michigan
% report bugs to:  hernan@umich.edu

mygray = [1:256]' * [1 1 1];

myhot = [128.5:0.5:256]' * [1 0 0] ;
tmp =   [128.5:0.5:256]' * [0 1 0] ;
tmp(:,1) = 256;
myhot = [myhot; tmp];
tmp =   [128.5:0.5:256]' * [0 0 1];
tmp(:,1:2) = 256;
myhot =  [myhot;  tmp;];
myhot = myhot(1:3:end, :);


myblue = myhot;
myblue(:,1) = myhot(:,3);
myblue(:,3) = myhot(:,1);

mygreen = myhot;
mygreen(:,1) = 0;
mygreen(:,2) = 256;
mygreen(:,3) = 0;%mygray(:,1);

mymap = [mygray; myhot ; myblue ; mygreen]/256;
colormap(mymap)

%%8/4/2014%%
%{
if computer=='PCWIN64'
    mymap = mymap(1:4:end,:);
    colormap(mymap);
end
%}
%%%

axis off
fprintf('\n ... color map reconfigured for overlays');
