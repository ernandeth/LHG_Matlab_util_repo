fname = 'vol0013.nii';

[d h] = read_nii_img(fname);
fd = zeros(size(d));

figure(1)
subplot(211), lightbox(d);
dfix = zeros(size(d));

for z=1:h.dim(4)
    fd = fft2(d(:,:,z));
    
    kphs = angle(fd);
    kplane = log(abs(fd));
    
    
    [xx yy] = meshgrid(linspace(-1,1,256), linspace(-1,1,256));
    
    x0 = ones(size(xx(:)));
    x1 = xx(:);
    x2 = xx(:).^2;
    x3 = xx(:).^3;
    x4 = xx(:).^4;
    
    x5 = xx(:).^5;
    x6 = xx(:).^6;
    x7 = xx(:).^7;
    
    y0 = ones(size(yy(:)));
    y1 = yy(:);
    y2 = yy(:).^2;
    y3 = yy(:).^3;
    y4 = yy(:).^4;
    
    y5 = yy(:).^5;
    y6 = yy(:).^6;
    y7 = yy(:).^7;
    
    P = [x0 x1 x2 x3 x4 x5 x6 x7   y0 y1 y2 y3 y4 y5 y6 y7 ];
    
    betas = pinv(P) * kplane(:);
    khat = P*betas;
    
    %khat = exp(khat);
    %kplane = exp(kplane);
    
    khat = reshape(khat,256,256);
    kres = (kplane - khat).^2;
    
    thresh= std(kres(:))*2.5;
    
    kfix = kplane;
    kfix(kres>thresh) = khat(kres>thresh);
    
    % fix the corners...
    kfix(1:20, 1:20) = kplane(1:20, 1:20);
    kfix(end-20:end, 1:20) = kplane(end-20:end, 1:20);
    kfix(1:20, end-20:end) = kplane(1:20, end-20:end);
    kfix(end-20:end, end-20:end) = kplane(end-20:end, end-20:end);
    
    figure(2)
    subplot(221), imagesc(khat), colorbar, title('poly fit')
    subplot(222), imagesc(kplane), colorbar, title('original kspace')
    subplot(223), imagesc(kres), colorbar, title('residuals')
    subplot(224), imagesc(kfix), colorbar, title('corrected')
    
    tmp = exp(kfix);
    tmp = tmp.* exp(i*kphs);
    tmp = abs(ifft2(tmp));
    dfix(:,:,z) = tmp;
    %
    
end

figure(1)
subplot(212), lightbox(dfix);

write_nii(['d_' fname], dfix, h, 0);
