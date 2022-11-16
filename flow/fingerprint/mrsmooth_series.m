function mrsmooth_series(nm);
% function mrsmooth_series(nm);
fprintf('\n... applying adaptive (markov filter) ...');

[raw hdr] = read_nii_img(nm);

for n=1:size(raw,1)
    tmp = reshape(raw(n,:), hdr.dim(2), hdr.dim(3), hdr.dim(4));
    tmp  = mrfilter(tmp, [3.75 3.75 6]);
    out(n,:) = tmp(:);
end

write_nii(['s_' nm], out, hdr, 0);

fprintf('... done');
return


