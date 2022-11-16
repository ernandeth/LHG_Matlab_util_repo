function mrsmooth_series(nm);

[raw hdr] = read_nii_img(nm);

for n=1:size(raw,1)
    tmp = reshape(raw(n,:), hdr.dim(2), hdr.dim(3), hdr.dim(4));
    tmp  = mrfilter(tmp, [10 10 10]);
    out(n,:) = tmp(:);
end

write_nii(['s_' nm], out, hdr, 0);

return


