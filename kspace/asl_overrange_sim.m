
for c=1:-0.01:0.1
    con = 100*phantom(128);
    tag = 0.99*con;
    
    nlevel = 10;
    kcon = fft2(con) + nlevel* randn(size(con));
    ktag = fft2(tag) + nlevel* randn(size(tag));
    
    con2 = ifft2(kcon);
    tag2 = ifft2(ktag);
    
    sub = con2-tag2;
    
    clipLevel = c *max(kcon(:));
    
    rktag = real(ktag);
    iktag = imag(ktag);
    
    rkcon = real(kcon);
    ikcon = imag(kcon);
    
    rktag(rktag > clipLevel) = clipLevel;
    rkcon(rkcon > clipLevel) = clipLevel;
    iktag(iktag > clipLevel) = clipLevel;
    ikcon(ikcon > clipLevel) = clipLevel;
    
    kcon = complex(rkcon, ikcon);
    ktag = complex(rktag, iktag);
    
    con3 = ifft2(kcon);
    tag3 = ifft2(ktag);
    
    sub3 = con3-tag3;
    
    subplot(221)
    imagesc(abs(con));
    subplot(222)
    imagesc(abs(con2));
    
    subplot(223)
    imagesc(angle(sub));
    %caxis([0 0.5])
    
    subplot(224)
    imagesc(angle(sub3));
    title(['clipping level : ' num2str(clipLevel)])
    %caxis([0 0.5])
    
    colormap gray
    drawnow
    
end
