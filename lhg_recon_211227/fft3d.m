function im = fft3d(kdata);
% function im = fft3d(kdata);
% does the fft in 3D and takes care of the fftshifts
% useful for 3D reconstruction

im=zeros(size(kdata));
xdim = size(kdata,1);
ydim = size(kdata,2);
zdim = size(kdata,3);

for sl=1:zdim
    im(:,:,sl) = fftshift(fft2(fftshift(kdata(:,:,sl))));
end
for x=1:xdim
    for y=1:ydim
        im(x,y,:) = fftshift(fft(fftshift(im(x,y,:))));
%        im(x,y,:) = (fft(fftshift(im(x,y,:))));
    end
end

return 