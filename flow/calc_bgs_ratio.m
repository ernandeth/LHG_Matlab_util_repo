[raw h] = read_nii_img('timeseries_mag.nii');

sup = 1e4 * raw(end,:) ./ raw(1,:);
h.dim(5) = 1;

write_nii('bgsratio.nii', sup, h,0);55