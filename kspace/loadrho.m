function rf = loadrho(fname)
% read .rho RF files in /usr/g/bin/ on scanner
% see also John Pauly's loadwave.m

fid = fopen(fname,'r','ieee-be');

% read entire file as int16
d = fread(fid,'int16');

% get RF waveform (real) 
% from John Pauly's loadwave.m, looks like last 4 shorts are part of header
rf = d(33:length(d)-4);

fclose(fid);

% display FT of pulse
m = fftshift(ifft(ifftshift(rf)));
%m = ifftshift(ifft(fftshift(rf)));

subplot(3,1,1);
plot(rf);
title('Magnitude')

subplot(3,1,2);
plot(abs(m));
title('FFT Magnitude')

subplot(3,1,3);
plot(unwrap(angle(m)));
title('FFT Phase (deg)')
%plot(phase(m));

% theta = unwrap(angle(m));
% %theta = theta(240:260);
% r = 152:168;
% theta = theta(r);
% 
% x = (r - r(1))';
% p = polyfit(x,theta,1)
% thfit = p(2) + p(1)*x ;
% figure;
% plot(theta-thfit);
% hold on;
%plot(thfit,'r');



return;
