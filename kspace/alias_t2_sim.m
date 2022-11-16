Npulses = 20;

myalphas = 10:5:30;

data = zeros(1000,2)
signal = zeros(Npulses,length(myalphas));

R1 = 1/1.2;
dt = 0.001;
a = 1;
TR = 30*dt;
E1 = exp(-TR*R1);



return
%% Aliasing
% This next simulation show swhat happens if the inversion plane starts
% wrapping around and getting into the image (Z-diraction aliasing)
% blurring effects are included too.
for pos=4:64
    
    % make an object
    im = zeros(128);
    im(33:96, 33:96) = phantom(64);
    % labeling band
    im(pos-2:pos+2,52:78) = 2;
    if 0  % this code changes the shape of the labeling band
        im(20:21,58:72) = 1;
        im(19,60:70) = 1;
        im(25:28,58:72) = 0;
        im(29,60:70) = 0;
    end
    
    % k-space:
    kim = (fft2(im));
    
    figure(1)
    subplot(221), imagesc(im); title('original')
    subplot(222); imagesc(abs(kim)); title('original kspace')
    
    % now undersample
    ukim = kim(1:2:end, 1:2:end);
    uim = fftshift(ifft2(ukim));
    
    subplot(223), imagesc(abs(uim)); title('undersampled');
    subplot(224); imagesc(abs(ukim)); title('undersmpled kspace')
    
    % now undersample with T2 weighting
    t2w = exp(-[1:128]* 0.19);
    t2w(65:end) = t2w(64:-1:1);
    t2w = kron(ones(1,128), t2w');
    
    wkim = t2w.*kim;
    
    uwkim = wkim(1:2:end, 1:2:end);
    uwim = fftshift(ifft2(uwkim));
    
    figure(2)
    subplot(221); imagesc(t2w); title('T2 weighting function')
    subplot(223), imagesc(abs(uwim)); title('undersampled and T2');
    subplot(224); imagesc(abs(uwkim)); title('undersmpled kspace and T2')
    subplot(222); plot(abs(uwim(:,32)),'r'); hold on; plot(abs(uim(:,32))); hold off; axis tight; title('vertical profil through the middle')
    
    drawnow
end


