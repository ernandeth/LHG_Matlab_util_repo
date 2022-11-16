function voxels = mask2voxel(mask)
% function voxels = mask2voxel(mask)
% convert from 3-D mask to voxel list in canonical orientation
% [i j k] = row, column, slice
% [x y z] in brain if brain is in analyze format
% (x is rows, y is columns, z is slices)
%
% Tor Wager, 10/17/01

voxels = [];
index = 1;
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        for k = 1:size(mask,3)
            if mask(i,j,k) > 0, voxels(index,:) = [i j k];,index = index+1;,end
        end
    end
end

return