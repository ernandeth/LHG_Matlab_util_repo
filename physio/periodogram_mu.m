function [P_y,f,Delta_f]=periodogram_mu(y,fs,freq_range)
% USAGE: [P_y,f,Delta_f]=periodogram_mu(y,fs,freq_range)
% INPUT: y is the ECG signal,
%        fs is the sampling frequency.
% OUTPUT: P_y is the periodogram computed over freq_range
%         which is default [0,60] Hz, fs is the sampling frequency.
% Magnus Orn Ulfarsson, 2007.
% Last modified: 01-31-07
if(nargin==2) 
    freq_range=[0;60];
end
T=length(y);
y=y-mean(y);
Fy=fft(y)/sqrt(T);
f=0:fs/T:floor(T/2)*fs/T;
[g,n_max]=min(abs(f-freq_range(2)));
f=f(1:n_max);
Fy=Fy(1:n_max);
P_y=abs(Fy).^2;
Delta_f=fs/2/T;
%-------------------------------------------------------