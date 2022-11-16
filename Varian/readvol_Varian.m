function myvol = readvol_Varian (dim)
% function myvol = readvol_Varian ([xdim ydim zdim] )

bits1 = 32;  % floating point data
endian = 'ieee-le';

myvol = zeros(dim(1), dim(2), dim(3));
slfiles = dir('slice*.fdf');

for sl=1:dim(3)
	slname = slfiles(sl).name;
	fid = fopen(slname , 'rb', 'ieee-le');
	disp('reading 2D...')
	status = fseek(fid, -dim(1)*dim(2)*bits1/8, 'eof');

	myvol(:,:,sl)=fread(fid,[dim(1), dim(2)],'float32',endian);
end

lightbox(myvol);
save volume myvol