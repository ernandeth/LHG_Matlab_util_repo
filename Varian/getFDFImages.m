% script to load in a whopping amount of varian fdf images
clear; close all;
% set up pamameters
location = '/home/histo/vnmrsys/data/s_2013061901/histoEpip_01.img/';
fileNames = dir([location '*.fdf']);
img1 = zeros(64,64,1000);
img2 = img1;
im3  = img1;
k = 1;
for u = 3:3:numel(fileNames)
img1(:,:,k) = getFDF(location, fileNames(u-2).name); %Slice 1 
img2(:,:,k) = getFDF(location, fileNames(u-1).name);  %Slice 2
img3(:,:,k) = getFDF(location, fileNames(u).name);  %Slice 3
k = k+1;
end