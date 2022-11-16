% Junk from CompCor.m

% JUNK, but maybe useful for comparison to other methods
%
% Optional: Run RETROICOR
% physdata = convertEXphysio('061112lh_phys_03',0.025)
% PhysioMat = mkASLPhysioMat('physio.dat',0.025, 2, Nslices, 4, 3.5);
% rmReg('2Dout', PhysioMat);
% !mv residuals.nii retro.nii
%
% 1. Construct Design Matrix
% 
% Nframes = 60;
% Nslices = 12;
% TR = 4;
% 
% onsets = [30:60:Nframes*TR];
% duration = 30;
% 
% reg = zeros(Nframes*TR,1);
% for o=1: length(onsets)
%     reg(onsets(o):onsets(o)+duration) = 1;
% end
% reg = conv(reg, spm_hrf(1));
% ASLreg = reg(1:TR:Nframes*TR);
% SUBreg = reg(1:TR*2:Nframes*TR);
% 
% DM = [SUBreg-mean(SUBreg)   ones(Nframes/2,1)];
% 
% ASLmodulation = ones(Nframes,1);
% ASLmodulation(1:2:end) = -1;
% rawDM = [ASLmodulation.*ASLreg ASLmodulation ones(Nframes,1) ASLreg];
% 
% figure
% imagesc(rawDM);

% New_rawDM = rawDM(1:2:Nframes, 1:3:4);

% aslsub(rootname1, 1, 1, Nframes, 0, 1, 0);
% CompCor_aslsub('tnoiseROI_residuals', 1, 1, Nframes, 0, 1, 0);
% 
% spmJr('sub', DM, [1 0]);
% CompCor_spmJr('CompCor_sub', DM, [1 0]);
% 
% figure
% no_cor = lightbox('Zmap_0001.img');
% figure
% comp_cor = lightbox('CompCor_Zmap_0001.img');
% 
% N_no_cor = reshape(no_cor,[64*64*Nslices 1]);
% 
% N_comp_cor = reshape(comp_cor,[64*64*Nslices 1]);
% 
% figure
% plot(N_no_cor, N_comp_cor, '*');
% xlabel('No Correction');
% ylabel('CompCor');
% title('Z score of CompCor vs. No Correction');
