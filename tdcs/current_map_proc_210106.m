%% compute B vector field map from current
% Using FSL for coregistration
% uses a non-stop time series
% current is oscillating
% this version does not use mcflirt, but FLIRT (mcflirt fails to coregister things properly)
close all

% units:  want everything in Gauss, meters, amps and seconds
TR = 6*1.25 ; % (3/20/19)
TR = 2; % 4/12/19
TR = 4; % 4/15/19
TR = 4*1.25 ; % (4/23/19)
TR = 2; % 5/6/19

TE = 10e-3; % s
TE = 30e-1; % 4/12/19
TE = 25e-3; % 4/15/19
TE = 16e-3; % 4/23/19
TE = 18e-3; % 5/6 / 19
TE = 40e-3; % 3/24/21


GAMMA = 42.56 * 1e6 / 1e4;  % gyromagentic constant in Hz/Gauss
u0 = 4*pi*10^(-7);          % unit=H/m. Permeability (electromagnetism)
I_freq =0.02; % square wave (3/20)
I_freq =0.02; % sine wave (4/20)
I_freq =0.05; % sine wave (5/6/19)

% frequency peak in spectrum:
Nyquist = 1/2/TR
peakF = I_freq / Nyquist
% Henry = Tesla*meter^2 /ampere
% so u0 units are  in Tesla*meter/ampere
% conversion:   u0 * 1e4 (G/T)
% so now it's in :   Gauss*m^2 / amp

u0 = u0 * 1e4;

XRES = 64;
Cscale = (XRES/64)^2;

%% recon

vols_re=[];
vols_im=[];
vols_abs=[];

recon = 0;
isSpiral = 1;
doCompCor = 0;
decorrelateCompcor =1;

if recon
    
    % clean up
    !rm *vol* *.txt
    !rm -r xforms
    !rm *.nii *.img *.hdr
    
    %% recon images and write them out as NIFTI
    % also keep a .mat file
    all_ims=[];
    pfiles=dir('P*')
    
    for n=1:length(pfiles)
        % get a  NIFTI header from the first image:
        if n==1
            sprec1(pfiles(n).name,  'N', 'l','q', 'fy','t', 1, 'n', XRES);
            %sprec1(pfiles(n).name,  'RR', 'N', 'l', 'fy','t', 1, 'n', XRES);
            volniifiles = dir('vol_e*.nii');
            [d h] = read_nii_img(volniifiles(1).name);
        end
        
        % recon complex images and concatenate the results into a time series:
        sprec1(pfiles(n).name,  'com', 'l', 'fy', 'n', XRES,'C',Cscale);
        volfiles=dir('vol_e*.mat');
        load (volfiles(1).name);
        
        % demodulate (recon error?)
        %
        demod = ones(size(imar));
        demod(1:2:end,:,:) = -1;
        demod(:, 1:2:end,:) = -demod(:, 1:2:end,:);
        
        % collected sagittal images, but need axials
        %{
        demod = permute(demod,[1,3,2]);
        imar = permute(imar,[1,3,2]);
        tmp = h;
        tmp.dim(3) = h.dim(4);
        tmp.dim(4) = h.dim(3);
        tmp.pixdim(3) = h.pixdim(4);
        tmp.pixdim(4) = h.pixdim(3);
        h = tmp;
        %}
        
        
        
        for m=1:length(volfiles);
            load (volfiles(m).name);
            
            % collected sagittal images, but need axials
            % imar = permute(imar,[1,3,2]);
            % imar = smooth3(imar,'gaussian',[5 5 5], 2) ;
            
            
            
            all_ims = [all_ims ; imar(:)'];
            write_nii(sprintf('vol_real_%04d.nii', m) ,  real(imar(:)) .* demod(:), h, 0);
            write_nii(sprintf('vol_im_%04d.nii', m) , imag(imar(:))  .*demod(:) , h, 0);
            write_nii(sprintf('vol_abs_%04d.nii', m) , abs(imar(:)), h, 0);
            
        end
        
        TDIM = size(all_ims, 1);
        
        !rm vol_e*.mat
        % save cplx_imgs.mat all_ims h
    end
    all_ims = all_ims';
    !fslmerge -t volseries_abs.nii  vol_abs*.nii
    !fslmerge -t volseries_real.nii  vol_real*.nii
    !fslmerge -t volseries_im.nii  vol_im*.nii
    
  %%  
    
    %% realignment of the real and imaginagry parts
    %{
    !mcflirt -in vol_real.nii -o rvol_real.nii  -cost corratio -plots -mats -refvol 3 -sinc_final -stats -dof 9
    !mcflirt -in vol_im.nii -o rvol_im.nii  -cost corratio-plots -mats -refvol 3 -sinc_final -stats -dof 9
    
    R = [];
    for t=0:TDIM-1
        Rmat = load(sprintf('./rvol_im.nii.mat/MAT_%04d', t));
        R = [R; Rmat(3,1:3)];
    end
    %}
    %% realignment using the magnitude images
    % we determine parameters from the magnitude images and apply them to
    % the real and imaginary components separately.
    
    volfiles=dir('vol_abs*.nii');
    ref_name = 'vol_abs_0003.nii';
   % ref_name = 'vol_abs_0280.nii';  %%%%  only for run 1 on 5/3/21 %%%%%
    
    
    R=zeros(length(volfiles), 3);
    mvt=zeros(length(volfiles), 12);
    
    parfor pos=1:length(volfiles)
        % for pos=[9] % 1:2:length(volFiles)
        
        tgt_name   =    sprintf('vol_abs_%04d.nii', pos);
        tgt2_name =    sprintf('vol_real_%04d.nii', pos);
        tgt3_name =    sprintf('vol_im_%04d.nii', pos);
        
        fprintf('\n coregistering target = %s to ref = %s ...\n\n' , tgt_name , ref_name);
        
        optstr = '-cost corratio   '
        
        % flirt_str =   sprintf(' flirt %s -ref %s -in %s  -out r%s -cost mutualinfo -omat xform_%04d.txt', optstr, ref_name, tgt_name, tgt_name, pos)
        % flirt_str =   sprintf(' flirt %s -ref %s -in %s  -out r%s -cost normmi -dof 9 -omat xform_%04d.txt', optstr, ref_name, tgt_name, tgt_name, pos)
        flirt_str =   sprintf(' flirt %s -ref %s -in %s  -out r%s -cost corratio -dof 9 -omat xform_%04d.txt', optstr, ref_name, tgt_name, tgt_name, pos)
        %flirt_str1 = sprintf('flirt -ref %s -in %s  -out r%s  -init xform_%04d.txt  -applyxfm -applyisoxfm 4',  ref_name, tgt_name, tgt_name, pos)
        flirt_str2 = sprintf('flirt -ref %s -in %s  -out r%s  -init xform_%04d.txt  -applyxfm -applyisoxfm 4',  ref_name, tgt2_name, tgt2_name, pos)
        flirt_str3 = sprintf('flirt -ref %s -in %s  -out r%s  -init xform_%04d.txt  -applyxfm -applyisoxfm 4',  ref_name, tgt3_name, tgt3_name, pos)
        
        system(flirt_str);
        system(flirt_str2);
        system(flirt_str3);
        
        Rmat = load(sprintf('xform_%04d.txt', pos));
        
        R(pos,:) =  Rmat(3,1:3);  % keep the rotation for the z component
        tmp = Rmat(1:3, :);
        mvt(pos,:) =tmp(:)';
        
        
    end
    
    save Movement.mat mvt R
    
    fprintf('\n... done with coregistration\n');
    fprintf('\nmerging raw... \n');
    !fslmerge -t volseries_real.nii vol_im_*.nii
    !fslmerge -t volseries_im.nii vol_real_*.nii
    
    [re h] = read_nii_img('volseries_real.nii');
    [im h]= read_nii_img('volseries_im.nii');
    
    c = complex(re,im)';
    write_nii('volseries_abs.nii', abs(c(:)), h, 0);
    h.datatype = 16;
    h.bitpix = 32;
    write_nii('volseries_ang.nii', angle(c(:)), h, 0);
    
    fprintf('\nmerging realigned ... \n');
    !fslmerge -t rvolseries_abs.nii  rvol_abs_*.nii
    !fslmerge -t rvolseries_real.nii  rvol_real_*.nii
    !fslmerge -t rvolseries_im.nii  rvol_im_*.nii
    
    %
    fprintf('\nsmoothing realigned ... \n');
    !fslmaths rvolseries_abs.nii -s 3 rvolseries_abs.nii
    !fslmaths rvolseries_real.nii -s 3 rvolseries_real.nii
    !fslmaths rvolseries_im.nii -s 3 rvolseries_im.nii
    %}
    fprintf('\ncleaning up ... \n');
    !rm *vol_im* *vol_real* *vol_abs*
    !mkdir xforms
    !mv xform_*.txt xforms/
    
else
    fprintf('\nLoading movement parms from .mat file ... \n');
    load Movement.mat
end

%read in the realigned images and prep them for processing
fprintf('\nReading realigned images ...\n');

[re h] = read_nii_img('rvolseries_real.nii');
[im h]= read_nii_img('rvolseries_im.nii');
all_ims = complex(re, im)';

write_nii('rvolseries_abs.nii' , abs(all_ims(:)), h, 0);

h.datatype = 16;
h.bitpix = 32;
write_nii('rvolseries_ang.nii' , angle(all_ims(:)), h, 0);

%% the current waveform
TDIM = h.dim(5);
t = 0:TDIM-1;
current_fun = sin(t*I_freq*TR*2*pi)';

%%
% make  current and baseline regressors for each position
% first identify the beginning and the end of each head postion
% you can spot the position changes because of large difference in
% movement vectors

z=zeros(TDIM,1);

% make a set of markers to indicate when the subject changes position
% use the change in the movement vectors to identify jumps.
mm = abs(diff(mean(mvt,2)));
markers = find(mm > 0.15 );
search = [1 markers' TDIM];

mvt2=[];
currents = [];
fin=0;
for n=1:length(search)-1
    
    beg = search(n);
    % does each period overlap with the previous period?
    if (fin >= beg)
        beg=fin+1;
    end
    if ( beg-fin > 1)
        beg=fin+1;
    end
    
    fin = search(n+1);
    
    % The different head positions are from beg to fin
    % [beg fin]
    
    rr = z;
    if (fin-beg) > 1  % make sure that they are not next door neightbors
        % make baseline regressor for each position:
        rr(beg:fin) = 1;  % DC component
        rr2 = rr;         % add a linear component
        rr2(beg:fin) = linspace(-1,1,fin-beg+1);
        mvt2 = [mvt2 rr rr2];
        %mvt2 = [mvt2 rr];
        
        % now make regressors for the current at each position
        cc = zeros(size(current_fun));
        cc(beg:fin) = current_fun(beg:fin);
        currents = [currents cc];
    end
    % add regressors to remove spikes at begining and end
    rr=z;
    rr(beg) = 1;
    mvt2 = [mvt2 rr];
    
    if (beg~=fin)
        rr=z;
        rr(fin) = 1;
        mvt2 = [mvt2 rr];
    end
    
end

% 2.14.2020 add the raw movement parameters to the postion markers
mvt3 = mvt;
for n=1:size(mvt,2)
    mvt3(:,n) =  mvt3(:,n) - mean(mvt3(:,n));
    mvt3(:,n) = mvt3(:,n)/max(abs(mvt3(:,n)));
end
mvt2 = [mvt2 mvt3];
% end 2.14.2020

imagesc(mvt2)


% 19.12.13:  include compcor clean up here
%
if decorrelateCompcor
    D=[currents mvt2];
    D = currents;
   
else
    D=[];
end

if doCompCor
    fprintf('\n Doing Compcor to clen up the realigned images ...\n');
    
    [re2 pts_junk] = compcor12(re, h, 20, D);
    [im2 pts_junk] = compcor12(im, h, 20, D);
    
    all_ims = complex(re2, im2)';
    h.datatype = 16;
    h.bitpix = 32;
    write_nii('rvolseries_ang_clean.nii' , angle(all_ims(:)), h, 0);
end
%
%% a model for the change in phase
% phase drift is modeled as linear

% the current function was a boxcar function:
%{
current_fun = ones(TDIM*TR*100,1);
period = 100/I_freq;
ons = period/2 : period : length(current_fun);
ons = ons; % + 2*TR*100 ;
dur=period/2;
for n=1:length(ons)
    current_fun( round(ons(n)) :  round(ons(n) + dur)) = -1;
end
current_fun = 0.5*smooth(current_fun,1000);
current_fun = current_fun(1:TR*100:end);
current_fun = current_fun(1:TDIM);
%}
%
% Current is a sine wave now instead:
t = 0:h.dim(5)-1;
current_fun = sin(t*I_freq*TR*2*pi)';

save current_fun.dat current_fun -ascii

N_original = size(re,1);

TDIM = size(all_ims, 2);
bad_images = [];
%{
%write_nii('rvol_abs.nii', abs(all_ims(:)), h, 0);
% fprintf('\nRemoving Corrupt images ...\n');
% figure out where the movement was too bad.
% we'll exclude these points from the model
ref_img = 3;
error = all_ims - repmat(all_ims(:,ref_img), 1, size(all_ims,2));
error = sum(abs(error) , 1);
% average RMS error per image
me = mean(error);
se = std(error);
% the bad images happen in the movement jumps
% look at the derivative (diff) of the error in time
ade = abs(diff(error));
sade = std(ade);
% create a vector indicating which images are the worst
badinds = zeros(size(error));
badinds(ade > sade/3) = 1;   % where the change in error is too high
bad_images = find(badinds);

% get rid of the neighbors, too
bad_images = [bad_images, bad_images+1 , bad_images-1];
bad_images(bad_images<1) = [];
bad_images(bad_images > length(badinds)) = [];
bad_images = sort(unique(bad_images));

% Not useful.... let's just not do this.
bad_images=[];


hold off
plot(error); hold on
plot(ade);
stem(badinds);
title('RMS difference in mag time series - image clean up')

% now remove the bad images (too much motion) from the series
all_ims(:,bad_images) = [];
R(bad_images,:) = [];   % remove the rotations from the bad images


%all_ims(:,bad_images) = 0;
%R(bad_images,:) = 0;   % remove the rotations from the bad images



% write out the new 'cleaned up' time series:
h.dim(5) = size(all_ims,2);
write_nii('rvolseries_abs_clean.nii', abs(all_ims(:)), h, 0);
h.datatype = 16;
h.bitpix = 32;
write_nii('rvolseries_ang_clean.nii', angle(all_ims(:)), h, 0);
fprintf('\n... Done Removing Corrupt images.\n');
%}
%%
fprintf('\nPhase unwrapping ....\n');
% temporal unwrapping
%
%ref_img = mean(all_ims(:,5:8),2);
ref_img = mean(all_ims(:,51:61),2);
pts = all_ims ./ repmat(ref_img, 1, size(all_ims,2));


pts = angle(pts);
pts = unwrap(pts,[],2);
pts(isnan(pts)) = 0;
pts(isinf(pts)) = 0;

h.datatype = 16;
h.bitpix = 32;
h.dim(5) = size(all_ims,2);
write_nii('rvolseries_ang_clean_unwrapped.nii', pts(:), h, 0);

%% do a filter;
w = [1 1 1]/3;
F = eye(TDIM);
for n=2:TDIM-1;
    F(n, n-1:n+1) = w;
end
pts = (F * pts')';

%% making a mask
msk = mean(abs(all_ims),2);
msk(abs(msk)< 0.25*mean(msk(:))) = 0;
msk(abs(msk)> 0) = 1;
msk(isnan(msk)) = 0;

% just the bottom slices for extracting the respiration
%{
msk2 = (reshape(msk, h.dim(2), h.dim(3), h.dim(4)));
msk2(:,:,11:end) = 0;
msk2 = msk2(:);
%}

figure; lightbox(reshape(msk, h.dim(2), h.dim(3), h.dim(4)));
for n=1:TDIM
    pts(:,n) = pts(:,n) .* msk;
end

Npositions = size(currents,2)

% X = [current_fun ones(TDIM, 1) linspace(-1, 1, TDIM)' mvt2];
% osc = ones(TDIM, 1); osc(2:2:end) = -1;

% X = [current_fun mvt2 ];
% currents = current_fun;  % disregard the splitting
X = [currents mvt2];

% mean center and scale the regressors
for r=1:size(X,2)
    tmp = X(:,r);
    %tmp = tmp - mean(tmp);
    tmp = tmp / max(abs(tmp));
    X(:,r) = tmp;
end

X(bad_images,:) = [];  % remove the parts of the model corresponding to the bad images
t(bad_images)= [];
X(isnan(X)) = 0;


figure;imagesc(X); title('Model')
save Model.mat X
imagesc(X)


% now estimate the model at each voxel
fprintf('\nEstimating each component of the model ....\n');
%invX = pinv(X);

% 1/6/21:  alternative is to use regularize least squares:
lambda = 500;
invX = inv(X'*X + lambda*eye(size(X,2)))*X';

betahats = invX * pts';
residuals = var(pts' - X*betahats);
c=eye(size(betahats,1));
% append another contrast at the top for the sum
% of the first 5 positions
c = [zeros(1,size(c,2)) ; c];
c(1,1:5) = 1;

% sanity check:  reconstruct time series from its components:
%yhat = X(:,1:Npositions)*betahats(1:Npositions,:);
yhat = X(:,1:5)*betahats(1:5,:);
yhat = yhat';
write_nii('yhat.nii', yhat(:), h, 0);

% calculate T statistics for each contrast
Tscore = zeros(size(c,1), size(betahats,2));

for n=1:size(betahats,1)
    for p=1:size(betahats,2)
        Tscore(n,p) = betahats(n,p) ./ sqrt(c(n,:)*invX*residuals(p)*invX'*c(n,:)');
    end
end

%% calculate an F statistic for the sum of currents
X_f = X;
X_r = X;
X_r(:,1:5) = [];

y_reduced   = pts' - X_r * pinv(X_r) * pts';
y_full      = pts' - X_f * pinv(X_f) * pts';
y_raw       = pts';

SSE_reduced = sum(y_reduced.^2,1);
SSE_full    = sum(y_full.^2,1);
SSE         = sum(y_raw.^2,1);

mn_reduced = mean(y_reduced,1);
mn_full = mean(y_full,1);

dof_full =      size(y_raw,1) - size(X_f,2);
dof_reduced =   size(y_raw,1) - size(X_r,2);
dof_SSE =       size(y_raw,1) - 1;

F = (SSE_reduced - SSE_full )/5 ./ (SSE_full /dof_full) ;
F_pval = 1-fcdf(F, 5, dof_full);


h.datatype = 16;
h.bitpix = 32;
h.dim(5) = size(betahats,1);
% save the coefficients for the different components
write_nii('betahats.nii', betahats, h, 0);

h.dim(5) = size(Tscore,1);
write_nii('Tscores.nii', Tscore, h, 0);
h.dim(5) = 1;
write_nii('residuals.nii', residuals, h, 0);

write_nii('Fscores.nii', F, h, 0);
write_nii('log10p.nii', -log10(F_pval), h, 0);

%% Mask to identify statistically significant voxels s
stat_msk = F_pval;
stat_msk(F_pval < 0.05) = 1;  % exclude voxels with too much noise/residual
stat_msk(F_pval >= 0.05) = 0;  % exclude voxels with too much noise/residual

stat_msk=1;

% clean up the data:
% now we need to remove linear , DC and movement terms:
fprintf('\nRemoving nuisance components from the time series ....\n');
Q = X(:, Npositions+1:end);

pts_clean_reg = (pts' - Q*pinv(Q)*pts');

%  remove the global mean signal from what's left
%{
fprintf('\nRemoving global mean from the time series ....\n');
gm = zeros(TDIM,1);
for n=1:TDIM
    gm(n) = mean(pts_clean_reg(n,:) .*msk2');
end

gm = gm -mean(gm);

pts_clean_reg = pts_clean_reg - gm*pinv(gm)*pts_clean_reg;
%}

%{ Excluding pixels with high residuals
%{
fprintf('\nExcluding pixels with high residuals ....\n');
msk3 = msk';
msk3(residuals > (mean(residuals) + 1.5*std(residuals))) = 0;
figure; lightbox(reshape(msk3, h.dim(2), h.dim(3), h.dim(4)));
for n=1:TDIM
    pts_clean_reg(n,:) = pts_clean_reg(n,:) .*msk3;
end
%}



h.dim(5) = size(pts,2);
tmp = pts_clean_reg';
write_nii('rvolseries_ang_regressed.nii', tmp(:), h, 0);


%{
% We can also try to reconstitute the data from the parameter estimates of
% the linear model
fprintf('\nReconstructing the time series from sine components....\n');
tmp = X(:,1:Npositions)*betahats(1:Npositions,:);
tmp = tmp';
write_nii('rvolseries_ang_compressed.nii', tmp(:), h, 0);
pts_clean_reg_filt = tmp';
%}

%% Calculate the delay in the waveform
%{
% compute the delay and amplitude of the wave from the first ten points of
% the time series.
%
% trig identity:
%
%   Bz*sin(w*t + phi) = Bz*sin(w*t) * cos(phi) + Bz * cos(w*t) + sin(phi)
% we have estimates for
%       beta1 = Bz*sin(phi)
%       beta2 = Bz*cos(phi)
% so ....
%       beta2/beta1 = tan(phi)
% and  ..
%       phi = atan2(beta1, beta2)
%       Bz = beta1/cos(phi)
%

% estimate the sine and cosine components from only the still period
% choose the first 16 scans
% PHI_pts = pts_clean_reg;
% PHI_X = X
% PHI_betahats = pinv(PHI_X) * PHI_pts';
% PHI_sigma2 = var(PHI_pts' - PHI_X*PHI_betahats);
% PHI_sigma2(PHI_sigma2>0.5) = nan;
% h2 = h; h2.dim(5) = 1;



% phi = atan2(PHI_betahats(1,:) , PHI_betahats(2,:) );
% % get rid of the symmetric wrap
% for n=1:length(phi)
%     if phi(n) < 0
%         phi(n) = pi-abs(phi(n));
%     end
% end

% find out which pixels had a signficant sinusoid
%Tscore =( PHI_betahats(1,:).^2 + PHI_betahats(2,:).^2 ) ./ PHI_sigma2;
%
%Tscore = PHI_betahats(1,:) ./ sqrt(PHI_sigma2);
%figure
%subplot(211);
%lightbox(reshape(Tscore, h.dim(2),h.dim(3),h.dim(4)),[-2 2], 20);
%caxis([-1 1])
%title('pixels with significant sinusoidal wave')

%
%
% phi = phi(Tscore>2);
%
% % figure out the delay from the most common phi
% [n val] = hist(phi,100);
% subplot(212); hist(phi,100);
% [mx bin] = max(n);
% imgPHI = val(bin);
% title(sprintf('phase estimates in sig. region: %f rads', imgPHI));

%mod_fun = sin(freq_n * t - imgPHI);
%Bz_amp= betahats(1,:) ./ cos(phi);
%}

%
% Test:  what happens if I don't regress out the movement - I think that it
% removes the projection info
% pts_clean = pts;

%% Now we solve for the  Bvector and the J vector
fprintf('\nSolving for the whole B field vector map ....\n');

mod_fun=2*current_fun';

% modulate the rotations by the sinusoid
Rmod = R .* repmat(mod_fun,3,1)';

%Bvec = pinv(Rmod) * pts_clean_reg;

% alternative:  use regularized least squares for the inversion:
Bvec = (Rmod'*Rmod + lambda*eye(size(Rmod,2)))*Rmod' * pts_clean_reg;

% from radians to Gauss
Bvec = Bvec/ (TE * GAMMA)  ;


%%
%3D smoothing - do it several times?
b_x = reshape(Bvec(1,:), h.dim(2), h.dim(3), h.dim(4));
b_y = reshape(Bvec(2,:), h.dim(2), h.dim(3), h.dim(4));
b_z = reshape(Bvec(3,:), h.dim(2), h.dim(3), h.dim(4));
%
% b_x = b_x .*msk;
% b_y = b_y .*msk;
% b_z = b_z .*msk;

for n=1:4
    b_x = smooth3(b_x,'gaussian',1) ;
    b_y = smooth3(b_y,'gaussian',1) ;
    b_z = smooth3(b_z,'gaussian',1) ;
end
%}


fprintf('\nSolving for the J field vector map from the curl of the B field....\n');

[xcoords , ycoords, zcoords] = meshgrid(...
    linspace(0, 1, h.dim(3)), ...
    linspace(0,1, h.dim(2)), ...
    linspace(0, 1, h.dim(4)) );

xcoords = xcoords * h.pixdim(2) *h.dim(2) * 1e-3;  % mm to meters
ycoords = ycoords * h.pixdim(3) *h.dim(3)  * 1e-3;  %
zcoords = zcoords * h.pixdim(4) *h.dim(4)  * 1e-3;  %


%
% calculate the current density from the magnetic field
% This will come out in units of Amps/m^2
% u0 is in Gauss*m2/amp

% note:   phantom:
% we used 2 mA through a tube of 0.25 cm radius
% 2e-3 / pi / 2.5e-3 / 2.5e-3
% so average should be ~ 10^2 A/m2

% human head:
%  2e-3 / pi / 0.07 / 0.07
% ~ ... ~ 0.13 A/m2

msk = reshape(msk,  h.dim(2), h.dim(3), h.dim(4));
[Jx, Jy, Jz] = curl(xcoords, ycoords, zcoords, b_x, b_y, b_z)  ;
Jx =  msk .* Jx / (u0);
Jy =  msk .* Jy / (u0);
Jz =  msk .* Jz / (u0);
Jmag =  (Jx.^2 + Jy .^2 + Jz.^2) .^0.5;

Jx = Jx(:).*stat_msk';
Jy = Jy(:).*stat_msk';
Jz = Jz(:).*stat_msk';
Jmag = Jmag(:).*stat_msk';



Bvec = Bvec .* repmat(stat_msk,3,1);

% write out components as separate files
h.dim(5) =1;

write_nii('Bx.nii',Bvec(1,:) , h, 0);
write_nii('By.nii',Bvec(2,:) , h, 0);
write_nii('Bz.nii',Bvec(3,:) , h, 0);
write_nii('B_mag.nii', (Bvec(1,:).^2 + Bvec(2,:) .^2 + Bvec(3,:).^2) .^0.5,  h ,0);

% write out components as separate files
write_nii('Jx.nii',Jx, h, 0);
write_nii('Jy.nii',Jy , h, 0);
write_nii('Jz.nii',Jz , h, 0);

write_nii( 'J_mag.nii', (Jx(:).^2 + Jy(:) .^2 + Jz(:).^2) .^0.5,  h ,0);
%%
%
ortho2005([],...
    'anat_file','volseries_abs.nii',...
    'wscale',[],...
    'spm_file','B_mag.nii',...
    'spm_file2',[],...
    'interact',1, ...
    'threshold',[1 1e3],...
    'threshold2',[],...
    'tseries_file','rvolseries_ang_regressed.nii',...
    'doFFT', 1 ...
);

%}rvolseries_ang_regressed