edit %   

% %  Field Maps of Shims

targetdir = '/home/histo/vnmrsys/data/s_2013011101';

cd(targetdir);

cd('gems_02.fid')

[imageMatrix1, kSpaceMatrixTmp1, procparTmp1] = reconstructImageArrayed('./');

% cd([targetdir '/gems_11.fid']);
% [imageMatrix2, kSpaceMatrixTmp2, procparTmp2] = reconstructImage('./');

angles1(:,:,1) = angle(imageMatrix1(:,:,1));
angles1(:,:,2) = angle(imageMatrix1(:,:,2));
angles1 = unwrap(angles1,[],3);

map1 = ((angles1(:,:,1)-angles1(:,:,2)));   % Units in gause
block1 = map1(26:41,26:39)/4257/0.002;
% block1 = unwrap(block1,[],2)/4257/0.002;
figure(1), imagesc(map1), title('Units in Gauss');  colorbar;
variance1 = std(block1(:));

%%
% angles2(:,:,1) = angle(imageMatrix1(:,:,2));
% angles2(:,:,2) = angle(imageMatrix2(:,:,2));
% angles2 = unwrap(angles2,[],3);
% 
% map2 = (angles2(:,:,1)-angles2(:,:,2))/4257/0.002;   % Units in gause
% block2 = map2(25:41,25:41);
% 
% figure(2), imagesc(block2), title('Units in Gauss'); colorbar;
% variance2 = std(block2(:));