%%% average matrices
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

figure, imagesc(N_av_c/5),colorbar, title('Necessity matrix (average)')
figure, imagesc(S_av_c/5),colorbar, title('Sufficiency matrix (average)')
figure, imagesc(rho_av_c/5),colorbar, title('Correlation matrix (average)')
