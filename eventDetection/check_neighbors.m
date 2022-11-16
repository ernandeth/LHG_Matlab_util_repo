function outdata = check_neighbors(indata, N, xdim,ydim,zdim)

outdata=zeros(size(indata));
Nframes=size(indata,1);
Npix=size(indata,2);
for pix=1:Npix
    [x y z] = ind2sub([xdim ydim zdim], pix);
    [nx, ny nz] = meshgrid([x-N:x+N], [y-N:y+N], [z-N:z+N]);
    
    % don't want stuff outside our volume
    inds = find(...
		nx(:)>0 & ...
		nx(:)<=xdim & ...
		ny(:)>0 & ...
		ny(:)<=ydim & ...
		nz(:)>0 & ...
		nz(:)<=zdim);
    % don't want the center pixel
    inds = inds(inds ~= floor(length(nx(:))/2)+1);

    nx2 = nx(inds);
    ny2 = ny(inds);
    nz2 = nz(inds);
    
    neighborinds = sub2ind([xdim,ydim,zdim], nx2(:), ny2(:), nz2(:));

	buffer = indata(:,neighborinds);
	outdata(:,pix) = sum(buffer,2);

end
    
outdata(outdata <= 1+N)=0;
outdata(outdata>0)=1;

outdata = outdata & indata;
            