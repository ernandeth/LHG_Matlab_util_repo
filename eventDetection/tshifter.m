function result = tshifter(y,tshift)
% function result = tshifter(y,tshift)
% 
% this function does a temporal shift of a signal by adding phase to it 
% in the frequency domain
Y = fftshift(fft(y));
Y = reshape(Y,length(Y),1);
w = linspace(-pi,pi,length(Y))';
Yshifted = exp(-i*w*tshift) .* Y ;
result = ifft(fftshift(Yshifted));
result = real(result);
return
