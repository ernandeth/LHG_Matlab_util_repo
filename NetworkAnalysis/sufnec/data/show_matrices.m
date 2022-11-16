i = 1; % subject 1->5

% load nocleanmat
load cleanmat
subplot(321),imagesc(Nclean{i}),colorbar, title('Nec (cleaned data)')
subplot(322),imagesc(Nnoclean{i}),colorbar, title('Nec')
subplot(323),imagesc(Sclean{i}),colorbar, title('Suf (cleaned data)')
subplot(324),imagesc(Snoclean{i}),colorbar, title('Suf')
load clean_rho
subplot(325),imagesc(rho{i}),colorbar, title('correlation (cleaned data)')
load noclean_rho
subplot(326),imagesc(rho{i}),colorbar, title('correlation')


%%% average matrices
figure
load cleanmat
N_av_c = 0;
for i=1:5
    N_av_c = N_av_c + Nclean{i};
end
S_av_c = 0;
for i=1:5
    S_av_c = S_av_c + Sclean{i};
end

load nocleanmat
N_av = 0;
for i=1:5
    N_av = N_av + Nnoclean{i};
end
S_av = 0;
for i=1:5
    S_av = S_av + Snoclean{i};
end


load clean_rho
rho_av_c = 0;
for i=1:5
    rho_av_c = rho_av_c + rho{i};
end
load noclean_rho
rho_av = 0;
for i=1:5
    rho_av = rho_av + rho{i};
end

subplot(321),imagesc(N_av_c/5),colorbar, title('Nec average (cleaned data)')
subplot(322),imagesc(N_av/5),colorbar, title('Nec average')
subplot(323),imagesc(S_av_c/5),colorbar, title('Suf average (cleaned data)')
subplot(324),imagesc(S_av/5),colorbar, title('Suf average')
subplot(325),imagesc(rho_av_c/5),colorbar, title('correlation average (cleaned data)')
subplot(326),imagesc(rho_av/5),colorbar, title('correlation average')
