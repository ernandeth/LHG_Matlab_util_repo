function nii2avw(name)


[dnii hnii] = read_nii_img(name);
avwh = nii2avw_hdr(hnii);

write_img([name '.img'],dnii,avwh);
return