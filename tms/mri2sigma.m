a = read_img2('t1spgr.img');
h = read_hdr('t1spgr.img');

% make isotropic matrix 64 x 64 x 64
[x,y,z] = meshgrid(1:h.xdim/64:h.xdim, 1:h.ydim/64:h.ydim, 1:h.zdim/64:h.zdim);

b = interp3(a,x,y,z);

h.xsize = h.xsize * h.xdim / 64;
h.ysize = h.ysize * h.xdim / 64;
h.zsize = h.zsize * h.xdim / 64;

h.xdim = 64;
h.ydim = 64;
h.zdim = 64;

h.origin = [32 32 32 0 0];

write_img('t1spgr64.img',b,h);
write_hdr('t1spgr64.hdr',h);

M = eye(4);
M(:,4) = [h.xsize*32; h.ysize*32; h.zsize*32; 1];
spm_segment('t1spgr64.hdr',eye(4),'');

% conductivities from Cerri et al 1996: J. Medical Engineering and
% Technology, 19, 1, p. 7-16
all = b;
% remove artifact:
all(:,3,:) = 0;
threshold = 0.5;
white = read_img2('t1spgr64_seg2.img');
white(find(white >= threshold)) = 1;
white(find(white < threshold)) = 0;
white = white * 1.2;

grey = read_img2('t1spgr64_seg1.img');
grey(find(grey >= threshold)) = 1;
grey(find(grey < threshold)) = 0;
grey = grey * 0.3;

csf = read_img2('t1spgr64_seg3.img');
csf(find(csf >= 0.2)) = 1;
csf(find(csf < 0.2)) = 0;
csf = csf * 1.6;

tmp = white + grey +csf;

bone = zeros(size(grey));
bone(find(all>=100)) = 1;
bone(find(tmp)) =0;
bone = bone * 0.01;

sigma = tmp + bone;

skin = 0.5;

save sigma.mat sigma