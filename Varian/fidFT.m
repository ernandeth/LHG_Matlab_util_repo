function im = fidFT(fname)
[pth , name, ext] = fileparts(fname);
procparname = fullfile(pth,'procpar');

[kspace,np,ntraces,nblocks] = ReadVarian2D(fname);
ns = ReadProcpar('ns',procparname);

k = reshape(kspace, ntraces, np/2);
k2  = reshape(k', ntraces/ns, np/2, ns);

im = zeros(size(k2));

for slice = 1:ns
    im(:,:,slice) = fftshift(fft2(fftshift(k2(:,:,slice))));
end
%{
subplot(121)
lightbox(angle(im));
subplot(122)
lightbox(abs(im));
drawnow
%}
save reconData
return
