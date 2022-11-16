function noise = makeAR1(len, rho)
% function noise = makeAR1(length, rho)
%
% returns a noise vector of 'length' smaples 
% with AR(1) characteristics, sampled at 'sampling_rate' (in Hz)
% the vector is normalized such that the variance is one and the 
% mean is zero



% 
% length = 1200;
% rho=0.9

noise = zeros(len,1);
noise(1)=randn(1,1);
for count=2:len
   noise(count) = noise(count-1)*rho + randn(1,1);
end

noise = noise/count;
% make sure that the variance is one and the mean is zero
% noise = noise - mean(noise);
noise = noise / sqrt(var(noise));

% subplot(211),plot(noise)
% subplot(212), plot(abs(fft(xcorr(noise))))
        
%         VarR = 0;
%         Vo = WKfun('AR+WN',rho,VarR,len);
%         % This is just a safety measure,to force homogeneous variance:
%         Vo = WKfun('Cov2Cor',Vo);
%         % Variance-covariance matrix
%         %V = Vo * noise_amp^2;
%         
%         noise = Vo * rand(len,1);;
return


