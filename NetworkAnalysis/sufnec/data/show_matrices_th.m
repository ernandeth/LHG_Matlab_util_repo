clear all
close all
clc

th =0;% 0.6;

load cleanmat
N_av_c = 0;
for i=1:5
    N_av_c = N_av_c + Nclean{i};
end
S_av_c = 0;
for i=1:5
    S_av_c = S_av_c + Sclean{i};
end
load clean_rho
rho_av_c = 0;
for i=1:5
    rho_av_c = rho_av_c + rho{i};
end

N = N_av_c/5;
S = S_av_c/5;
R = rho_av_c/5;

Nbin = N>=th;
Sbin = N>=th;
N = N.*Nbin;
S = S.*Sbin;

%{
subplot(131),imagesc(N),colorbar, title('Necessity average'), axis image
subplot(132),imagesc(S),colorbar, title('Sufficiency average'), axis image
subplot(133),imagesc(R),colorbar, title('Correlation average'), axis image
%}

figure(1),imagesc(N),colorbar, title('Necessity average'), axis image
figure(2),imagesc(S),colorbar, title('Sufficiency average'), axis image
figure(3),imagesc(R),colorbar, title('Correlation average'), axis image
