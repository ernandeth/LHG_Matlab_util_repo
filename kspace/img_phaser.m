function img_phaser(root, ph)
% function img_phaser(root, ph)

names = dir([root '*.mat']);
Nframes = length(names);

load(names(1).name);

h = define_avwhdr;
h.xdim = size(imar,1);
h.ydim = size(imar,2);
h.zdim = size(imar,3);
h.tdim = Nframes;
h.xsize = 3.75;
h.ysize = 3.75;
h.zsize = 7;
h.datatype = 4;
h.bits = 16;

imar = imar(:);
data = zeros(Nframes, length(imar));

for t=1:Nframes
	load(names(t).name);
	imar = imar(:);
	data(t,:) = imar * exp(i*ph);
end


write_img('phased.img', abs(data), h);
write_img('p_phased.img', 1000*angle(data), h);

fprintf('results are saved as phased.img and p_phased.img');

