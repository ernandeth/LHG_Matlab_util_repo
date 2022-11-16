% Junk from tnoiseROI.m

%   JUNK but maybe useful later
%
%
%  rootname2: RETROICOR processed ANALYZE file 
% PhysioMat = mkASLPhysioMat('physio.dat',0.025, 2, 12, 4, 3.5);
% load PhysioMat.mat                
% polyMat = PhysioMat(:,:,1:4);
% t = (0:Nframes-1)';
% polynomial = [ones(size(t)) t/sum(t) t.^2/sum(t.^2)  t.^3/sum(t.^3) ];
% 
% rmReg(rootname1, polynomial);
% 
% [Y1, h1] = read_img('raw_residuals.nii');
% [Y2, h2] = read_img(rootname2);
% NewY2 = Y1-Y2;

%    Fractional Variance of Physiological Noise
%    - The ratio of the variance of the voxel time series determined with
%    RETROICOR to the variance of the original time series after removal of
%    constant and linear trends (the square of tSTD)

%  5.Plot tSTD versus fractional variance
% figure
% plot(Frac_VAR,tSTD,'*');
% title('tSTD vs. Fractional Variance of Physiological Noise');
% VAR_RETRO = zeros(n,1);
% Frac_VAR = zeros(n,1);
% VAR_RETRO(k,1) = var(NewY2(1:m,k));
% Frac_VAR = VAR_RETRO./(tSTD.^2);

% I1 = find(isinf(Frac_VAR));
% I2 = find(isnan(Frac_VAR));
% Frac_VAR(I1,1) = 0;
% Frac_VAR(I2,1) = 0;

% [r, s] = max(tSTD);
% 
% t = 0:m-1;
% t = t*TR;
% hi_t = Y1(:,s);
% hi_f = abs(fftshift(fft(hi_t)));
% 
% omega = [-1/(2*TR):(1/TR)/m:(1/(2*TR)-(1/TR)/m)];
% 
% figure
% subplot(2,1,1)
% plot(t, hi_t);
% subplot(2,1,2)
% stem(omega,hi_f);


%  6.Plot Fractional Variance of Physiological Noise vs. % of Voxels

% [sort_Frac_VAR, I1 ]= sort(Frac_VAR,'descend');
% 
% %    a. Determine how many different values are in Frac_VAR
% x = 1;
% 
% for k = 2:n
%     if (sort_Frac_VAR(k,1) < sort_Frac_VAR(k-1,1))
%         x = x+1;
%     end
% end

% %   b. Create matrices for plot
% FVP = zeros(x,1);
% PV  = zeros(x,1);
% FVP(1,1) = sort_Frac_VAR(1,1);
% 
% y = 1; % Counter for how many values in Frac_VAR are the same
% z1 = 1; % Counter for FVP index
% z2 = 1; % Counter for PV index
% 
% for k = 1:n 
%    if (k ~= n) && (sort_Frac_VAR(k,1) == sort_Frac_VAR(k+1,1))
%        y = y+1;  
%    else
%        FVP(z1,1) = sort_Frac_VAR(k,1);
%        PV(z2,1) = (y/n)*100;
%        z1 = z1+1;
%        z2 = z2+1;
%        y = 1;
%    end
% end
% 
% per_vox = cumsum(PV);
% 
% figure
% plot(per_vox,FVP);
% title('Fractional Variance of Physiological Noise vs. % of voxels');
% grid;
% axis([0 20  0.4 max(FVP)]);


% Display voxels with fractional variance value less than 0.3 and tSTD
% above 150

% for k = 1:n
%     if (tSTD(k,1) < 150) || (Frac_VAR(k,1) >= 0.3)
%         NewY3(:,k) = 0;
%     end
% end
% 
% 
% X1 = reshape(NewY3(1,1:49152),[64,64,12]);
% X2 = reshape(NewY3(2,1:49152),[64,64,12]);
% X3 = reshape(NewY3(3,1:49152),[64,64,12]);
% X4 = reshape(NewY3(4,1:49152),[64,64,12]);
% X5 = reshape(NewY3(5,1:49152),[64,64,12]);
% X6 = reshape(NewY3(6,1:49152),[64,64,12]);
% 
% NewX1 = zeros(64,64);
% NewX2 = zeros(64,64);
% NewX3 = zeros(64,64);
% NewX4 = zeros(64,64);
% NewX5 = zeros(64,64);
% NewX6 = zeros(64,64);
% 
% for k = 1:12
%     NewX1 = X1(:,:,k)+NewX1; 
%     NewX2 = X2(:,:,k)+NewX2; 
%     NewX3 = X3(:,:,k)+NewX3; 
%     NewX4 = X4(:,:,k)+NewX4; 
%     NewX5 = X5(:,:,k)+NewX5; 
%     NewX6 = X6(:,:,k)+NewX6; 
% end
% 
% 
% figure
% 
% subplot(2,3,1)
% imagesc(NewX1);
% 
% subplot(2,3,2)
% imagesc(NewX2);
% 
% subplot(2,3,3)
% imagesc(NewX3);
% 
% subplot(2,3,4)
% imagesc(NewX4);
% 
% subplot(2,3,5)
% imagesc(NewX5);
% 
% subplot(2,3,6)
% imagesc(NewX6);

% display_slice = zeros(64,64,12,60); % image width, image height, slice, time
% 
% dim1 = [64:64:49152];
% dim2 = 1;
% 
% for k = 1:64
%     if k == 1
%         display_slice(k,:,1,1) = NewY1(1,1:dim1(k));
%         dim2 = dim1(k)+1;
%     else
%         display_slice(k,:,1,1) = NewY1(1,dim2:dim1(k));
%         dim2 = dim1(k)+1;
%     end
% end
% 
% figure
% imagesc(display_slice(:,:,1,1));


% stim_sig1 = rawDM(:,1);
% stim_sig2 = rawDM(:,2);
% stim_sig3 = rawDM(:,3);
% stim_sig4 = rawDM(:,4);
% 
% for k = 1:z3
%     [RHO PVAL] = corr(NewY1(:,ind_noise_vox(k,1)), stim_sig1);
%     if PVAL < 0.2 % exclude voxels with a p-value less than 0.2
%         NewY1(:,ind_noise_vox(k,1)) = 0;
%     end
%     [RHO PVAL] = corr(NewY1(:,ind_noise_vox(k,1)), stim_sig2);
%     if PVAL < 0.2 % exclude voxels with a p-value less than 0.2
%         NewY1(:,ind_noise_vox(k,1)) = 0;
%     end
%     [RHO PVAL] = corr(NewY1(:,ind_noise_vox(k,1)), stim_sig3);
%     if PVAL < 0.2 % exclude voxels with a p-value less than 0.2
%         NewY1(:,ind_noise_vox(k,1)) = 0;
%     end
%     [RHO PVAL] = corr(NewY1(:,ind_noise_vox(k,1)), stim_sig4);
%     if PVAL < 0.2 % exclude voxels with a p-value less than 0.2
%         NewY1(:,ind_noise_vox(k,1)) = 0;
%     end
% end