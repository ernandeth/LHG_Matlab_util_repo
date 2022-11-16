function im = fidFT3d(fname)
[pth , name, ext] = fileparts(fname);
procparname = fullfile(pth,'procpar');

[kspace,np,ntraces,nblocks] = ReadVarian2D(fname);
ns = ReadProcpar('nv2',procparname);
ydim = ReadProcpar('nv',procparname);
xdim = np/2;
im = zeros(xdim,ydim,ns);

k = reshape(kspace,ydim,ns, xdim);

im = fftshift(fftn(fftshift(k)));

subplot(121)
lightbox(angle(im));
subplot(122)
lightbox(abs(im));
drawnow

save reconData
return