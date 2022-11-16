function FF=FourierMat(N)
% function F=FourierMat(N)
%
% builds a complex Fourier coefficient matrix of size N
% you can do and FFT of x by doing:  F*x
%

w = linspace(0, 2*pi - 2*pi/N, N);
F = zeros(N);
t = [0:N-1];
for n=1:N
    F(n,:) = exp(-i*w(n)*t); 
end

% turns out matlab already had that built in!
% F = dftmtx(N);

FF=zeros(N*N);
for n=1:N
    for m=1:N
        FF( (n-1)*N+1: n*N, (m-1)*N+1: m*N) = F * exp(-i*w(n)*t(m));
    end
end
% 
% 
% imagesc(angle(FF))
% 
% x=phantom(N);
% subplot(221), imagesc(x);
% x= x(:);
% xf = FF*x;
% subplot(222), imagesc(abs(reshape(xf,N,N)))
% xff = FF*xf;
% subplot(223), imagesc(abs(reshape(xff,N,N)))



return