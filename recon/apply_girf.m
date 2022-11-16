function corrected = apply_girf(mygrad, realp, imagp, cutoff, dt)
% function corrected = apply_girf(mygrad, realp, imagp, dt)
%
% mygrad    nominal input gradient
% realp  polynomial coefficients for real part of teh transfer frunction
% imagp  polynomial coefficients for imaginary part of teh transfer frunction
% cutoff maximum frequency in the transfer function
% dt     sampling interval of the nominal input gradient
%

% Define the frequencies over which we compulte the transfer function
% (just the positive side of the freq. spectrum)
Nyq = 0.5/dt;
df = 1/(dt*length(mygrad));
Npoints = length(mygrad);
f = linspace(0, Nyq-df, Npoints/2);

% The coefficients to the real part of the transfer function:
p = realp;
% reconstruct the transfer function at this sampling rate:
realjunk = polyval(p,f);
% but the girf is known over a different frequency range, so we zero out
% the unknown stuff


sigma = cutoff/50;

Ncutoff = round(length(realjunk)*cutoff/Nyq);
realjunk(Ncutoff+1 : end) = 0;

% then we rebuild the negative part of the spectrum by assuming symmetry
% between positive and negative parts of the spectrum
% (remember:  I havent' done any fftshifts)
realjunk = [ realjunk  0 realjunk(end:-1:2)];

realjunk = gaussiansmooth(realjunk, sigma, 5);

% Now the same for the imaginary part
p = imagp;
imjunk = polyval(p,f);

imjunk(Ncutoff+1 : end) = 0;

imjunk = [imjunk 0 -imjunk(end:-1:2) ];

imjunk = gaussiansmooth(imjunk, sigma, 5);

% make the 'ideal' transfer function from its real and imaginary parts
H_nice = complex(realjunk, imjunk);

% % % %%%%%%%%% a delay fudge here %%%%%%%%%%%%%%%
% phi = linspace(-Nyq, Nyq-df, Npoints);
% mydel = exp(-i*phi*2*pi* 21e-6);
% H_nice = H_nice .* mydel;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upsample the Transfer function to the appropriate sampling rate:
zpad = zeros(1, Npoints-length(H_nice));
H_nice = [H_nice(1:end/2) zpad H_nice(end/2+1 : end) ];

% the gradient's impulse response function in the imte dpmain
% is the iFFT of the transfer function in frequency

girf = ifft(H_nice).';
girf = girf(1:Npoints/2-1);

%%%
girf = ifftshift(ifft(H_nice)).';


%%%

% finally, the corrected waveform is the convolution with the GIRF
corrected = conv(mygrad, (girf));
%corrected = real(corrected(1:length(mygrad)));
corrected = real(corrected(length(girf)/2 : end - length(girf)/2) );

% (show me)
figure(76)
subplot(221)
plot(ifftshift(abs(H_nice)));
title('Magn. of Transfer fucntion')
subplot(223)
plot(ifftshift(angle(H_nice)));
title('Phase of transfer function')

subplot(222)
plot(real(girf)); hold on
plot(imag(girf),'r')
hold off
title('input response function')
legend('real', 'imaginary')

subplot(224)
plot(mygrad); hold on
plot(corrected,'r'); hold off
title('Gradient Waveform')
legend('nominal', 'corrected')

return

function  out = gaussiansmooth(in, sigma, Npts)


gkernel = 1/sigma/sqrt(2*pi) * exp( - (linspace(-5*sigma, 5*sigma, Npts).^2 /(2*sigma^2) ) );
gklen = length(gkernel);

junklen = length(in);
junknrg = norm(in);

out = conv(in, gkernel);
out = out(gklen/2+1:end-gklen);
out = out * junknrg / norm(out);

return


%% some code to play with
load girf_parms
allrealp
allimagp
allH_nice
dt
Nyq
cutoff
allgirf

traj = linspace(0,1,5000).^0.8.*exp(i*linspace(0,50*pi,5000));
plot(traj)

mygrad = real(traj);
dt = 1e-5;
corrected = apply_girf(mygrad, allrealp(1,:), allimagp(1,:), cutoff, dt);
