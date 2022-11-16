function phase_shifter(root, phase)
%function phase_shifter(root, phase)

[raw, h] = read_img(root);
h = nii2avw_hdr(h);
praw = read_img(['p_' root]);

raw = raw.* exp(-i .* praw/1000);

% here is the unwrapping action : add some phase to everything!

raw = raw.* exp(-i .* praw/1000);

write_img('p_shifted.img',1000*angle(raw),h);
write_img('shifted.img',abs(raw),h);

return