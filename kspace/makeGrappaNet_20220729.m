function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
% function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
%
% this version tried to convert neighbor locations into weights.
% It grids the calibration data into a Cartesian space first
% then calculates grappa weights from those data
% then trains a network to solve for those weights , given locations of neighbors
%
% Hopefully the next version will not need gridding to cartesian space first?
%
% first determine a center region for calibration assuming that this
% center region is well sampled.
%
% kx, ky, kz is the trajectory (cm*-1)
% data is the signal
% dim3 is num. pixels in each dimension of IMAGE space
% CalRadius is the fraction of k space used for calibration region
%
% this works up to the point of estiamting weights frorm locations.  Do not
% touch!!!

R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
Ncoils = size(data,2);

% the calibration region is determined by CalRadius relative to the Rmax
dim3_cal = floor(dim3*CalRadius)+1;  % Now dim3 will be the num. pixels in cal region.
dim3_cal = 2*floor(dim3_cal/2) +1;   % Making sure it's an odd number

CalRadius = CalRadius * Kmax;  % convert from fraction to k space units.
CalRegion = find(R<CalRadius); % indices of data in the calibration region

%%
fprintf('\rGRIDDING: Making cartesian calibration data from non cartesian data...\n')
% now interpolate spiral data into a cartesian grid
% will substitute this with kbgrid eventually
[subkx , subky, subkz ]= meshgrid( ...
    linspace(-CalRadius,CalRadius,dim3_cal(1)), ...
    linspace(-CalRadius,CalRadius,dim3_cal(2)), ...
    linspace(-CalRadius,CalRadius,dim3_cal(3)));

% dimensions: (calregion_Npix x Ncoils)
cal_grid = [];
dcf = [];
im = zeros(size(subkx));
for coilnum=1:Ncoils
     
    % calibration region for each coil
    fprintf('\rcoil  %d ... ', coilnum);
    tic
    %
    [tmp dcf kernel] = grid3d_lhg(...
        kx, ky, kz, ...
        data(:,coilnum), ...
        subkx, subky, subkz, ...
        2, dcf);
    toc
    %}
    %tmp = griddata(kx, ky, kz, data(:,coilnum), subkx, subky, subkz,'linear');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    %{
    % testing images for recon
    im = im + (fftshift(abs(fftn(tmp)))).^2;
    fk = fftshift(abs(fftn(kernel)));
    
    subplot(3,1,1)
    lightbox(log(abs(tmp )),[],5,[]); 
    title(['coil ' num2str(coilnum) ' ACS region interpolated']); 
    
    subplot(3,1,2)
    lightbox(sqrt(im(10:end-10,10:end-10, 10:end-10)),[],5,[]); 
    %lightbox(sqrt(im),[],4,[]); 
    title(['coil ' num2str(coilnum) ' image space']); 
    
    subplot(3,1,3)
    lightbox(abs(fk(10:end-10,10:end-10, 10:end-10)),[],5,[]); 
    title('kernel intensity artifact');
    
    drawnow
    %}
    cal_grid = [cal_grid tmp(:)];
end
%{
% testing the reconned images
figure(37)
fk = fk- min(fk(:));
fk = fk/max(fk(:));
im_all = sqrt(im)./(0.01+fk);
lightbox(im_all(10:end-10, 10:end-10, 20:end-20));
title('Final Corrected Image');
%}

% size of k-space voxels: distance between neighbors in the iso-tropic cartesian grid
dkx = abs(subky(2,1,1) - subky(1,1,1));

save workspace.mat
%%
fprintf('\nConstructing patches within calibration region ...\n ');
load workspace.mat

doRandomNeighbors = 1;

Nnbrs = 21;     % number of k space data locations for interpolation
% in each patch.  a patch is a constellation of points in
% kspace arround a location that we want to interpolate.

pklocs = [];    % The kspace coordinates of the data in each patch
% dimensions are: 3 x (Nequations-1)

Npatches = 1000;  % number of times we'll calcuate the interpolation kernel



% vector with possible kspace targets (the middle half of the calibration
% region).  This is also the allowed kspace distance to possible neighbors
% in units of pixels (k-space)
%
center_coords = [-floor(dim3_cal(1)/16):floor(dim3_cal(1)/16)];
nbr_coords = [-floor(dim3_cal(1)*7/16):-1   1:floor(dim3_cal(1)*7/16)];


if doRandomNeighbors==0
    nbr_coords = [ -2 -1 1 2 ];
    Nnbrs = length(nbr_coords)^3;
end


%  Targets: coordinates for the center of the patches
[x y z] = meshgrid(center_coords, center_coords, center_coords);
x = x(:);
y = y(:);
z = z(:);

% Neighbors:  coordinates for the neighbors.  They are equally spaced
% on the neighborhood around the center
% (use the neighbors within the dim3_cal/8 region of the target)
[dx dy dz] = meshgrid(nbr_coords, nbr_coords, nbr_coords);
dx = dx(:);
dy = dy(:);
dz = dz(:);

%Nequations = round(Ncoils*Nnbrs*1.5);
Nequations = length(x);
% Number of equations in a patch.
% The number of unknowns in the interpolation kernel
% is Ncoils*Nnbrs,
% so we need AT LEAST Ncoils*Nbrs  equations
% to solve that many unknowns per grappa patch


% GRAPPA interpolation does this:
%
% Ktarget = Weights .* Kneighbors
%
% Kneighbors    % dimensions: Nequations x Ncoils*Nnbrs
% Weights       % dimensions: Ncoils*Nnbrs x Ncoils
% Ktargets      % dimensions: Nequations x Ncoils
%

Kneighbors = zeros(Nequations, Ncoils*Nnbrs);
Ktargets = zeros(Nequations, Ncoils);

TrainDataIn = zeros( Ncoils*Nnbrs , 3, 1, Npatches*Nequations);
% TrainDataOut = zeros( Nequations*Npatches, Ncoils*2);
% we're going to do only one coil at a time
TrainDataOut = zeros( Nequations*Npatches, Ncoils*Nnbrs*2 );

train_cnt = 1;

for p = 1:Npatches
    
    if doRandomNeighbors
        % the neigbors are randomized for this patch
        % k-neighborhood for interpolation
        % units of pixels (k-space)
        dx = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
        dy = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
        dz = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    end
    
    nbrlocs =  [dx dy dz];
    nbrlocs = repmat(nbrlocs, Ncoils,1);
    
    Kneighbors(:) = 0;
    % COnstruct system of equations to solve for GRAPPA weights
    for n=1:Nequations
        % choose target at random to serve as center of patch
        ind = randi(length(x), 1);
        
        patchcenter = [x(ind) y(ind) z(ind)];          % grid location of the target
        pklocs =  [x(ind)+dx  y(ind)+dy z(ind)+dz];    % grid locations of the neighbors
        
        % put all the coordintates into one matrix
        % the first k space location is the center of the patch
        pklocs = [patchcenter; pklocs];
        
        % the k space locations of interpolation neighbors are relative to
        % the center of the calibraiton region
        pklocs = pklocs + ceil(dim3_cal/2);
        
        % Extract the signals at those locations
        % (remember that the first one is the center/target)
        inds = sub2ind(dim3_cal, pklocs(:,1), pklocs(:,2), pklocs(:,3));
        
        tgt = cal_grid(inds(1),:);  % These are the target complex signals
        % dimensions: 1 x Ncoils,
        Ktargets(n,:) = tgt;        % ... we stack them into Nequations x Ncoils
        
        % Neighbor complex signals (Nnbrs x Ncoils)
        tmp = cal_grid(inds(2:end),:);
        Kneighbors(n,:) =  tmp(:)';
        
        % store these for training the neural net
        % make a table of locations and signals for this target equation.
        
        %       nbrs_info = [nbrlocs real(tmp(:)) imag(tmp(:)) ];
        %       TrainDataIn( :, :, 1, train_cnt) = nbrs_info;
        
        TrainDataIn( :, :, 1, train_cnt) = nbrlocs;
        
        train_cnt = train_cnt +1;
        
    end
    
    
    % Solve for the GRAPPA weights.
    fprintf('\rSolving GRAPPA weights for patch ... %d ', p);
    
    [W, rnk, RMSE] = calc_grappa_kernel(Ktargets, Kneighbors);
    % Solve a system of equations, in which the same set of weights 
    % is used to interpolate the whole set of targets, since the neighbors are 
    % at the same relative positions. (It's a convolution kernel)
    % 
    % dimensions:
    % Kneighbors :  Nequations x (Nnbrs*Ncoils)
    % W:            (Nnbrs*Ncoils) x Ncoils
    % Ktargets :    Nequations x Ncoils
    
    % store the error for coil 1
    rmse(p) = RMSE(1);
    
    % The output of the Neural net will be W : 
    % the GRAPPA weights for each neighbor and each coil 
    % needed to interpolate a single target (for each coil)
    % (Nnbrs*Ncoils) x Ncoils
    
    % HOWEVER ... We'll start trying to do just one coil at a time
    coilnum = 1;
    Wblock = repmat(W(:,coilnum)', Nequations,1);
    Wblock = [real(Wblock) imag(Wblock)];
    
    % weed out training data if the GRAPPA fit is not very good
    if mean(RMSE) > 1e-8
        Wblock(:) = nan;
        fprintf('RMSE =%f too high - discard', mean(RMSE))
    end
    
    TrainDataOut(Nequations*(p-1)+1:Nequations*p,:) =  Wblock;
    %{
    if mod(p,20)==0
        subplot(221)
        imagesc(abs(Kneighbors))
        colorbar
        title('neighbors')
        xlabel('N coils * N neighbors')
        ylabel('N targets')
        
        subplot(222)
        imagesc(abs(W))
        colorbar
        title('Weights')
        xlabel('N coils ')
        ylabel('N neighbors* N coils')

        subplot(223)
        imagesc(abs(Ktargets))
        colorbar
        title('targets')
        xlabel('N coils ')
        ylabel('N targets')
        
        
        subplot(224)
        plot(RMSE)
        title(sprintf('Error, rank(Nbrs): %d ', rnk))
        xlabel('coils')
        ylabel('NRMSE (%)')
        
        drawnow
    end
    %}
    
end


% remove Training data containing NaNs
[a,b] = find(isnan(TrainDataOut));
TrainDataIn(:,:,:,a) = [];
TrainDataOut(a,:) = [];


% split the training data into training and testing data sets.
Ntrain = size(TrainDataIn,4);
%
% Grab a fraction of them for testing
inds = randperm(Ntrain, ceil(Ntrain/20));
TestDataIn = TrainDataIn( :,:,1 , inds);
TestDataOut = TrainDataOut(inds, :);

% ... and remove those from the training
TrainDataIn( :,:,:, inds) = [];
TrainDataOut(inds, :) = [];
save TrainData.mat TrainDataIn TrainDataOut TestDataIn TestDataOut
%
% construct a neural net that will take Neighbor intensities and locations and return the
% corresponding interpolated targets
% input dimensions:  (Nnbrs * Ncoils x 5)
% output dimensions:  1x2 ---> will be Ncoils x 2
%%
load TrainData.mat
GWnet = [
    
imageInputLayer([Ncoils*Nnbrs 3]);
%convolution2dLayer(3,5,'Padding','same');
%fullyConnectedLayer(Ncoils);

fullyConnectedLayer(Ncoils*Nnbrs*6);
batchNormalizationLayer;
reluLayer;
%tanhLayer;

% fullyConnectedLayer(Ncoils);
% batchNormalizationLayer;
% reluLayer;
%
% fullyConnectedLayer(Ncoils);
% batchNormalizationLayer;
% reluLayer;

% one coil only
% fullyConnectedLayer(2*Ncoils);
fullyConnectedLayer(2*Nnbrs*Ncoils);
%reluLayer
regressionLayer
]

scale = max((TrainDataOut(:)));
scale = 1;



TrainDataOut = TrainDataOut/scale;
TestDataOut = TestDataOut/scale;

fprintf('\nTraining Network ...\n ');
options = trainingOptions('sgdm',...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-4, ... %1e-4
    'WorkerLoad', ones([40 1]), ...
    'MaxEpochs',20, ...
    'LearnRateSchedule', 'piecewise', ...
    'Plots','training-progress');  % 'training-progress'

[GWcalcNet, info] = trainNetwork(TrainDataIn, TrainDataOut, GWnet, options);

save GRAPPAnet.mat GWcalcNet info scale
%%

fprintf('\nTesting the network with data from calibration region...')
justafew=1:100:size(TestDataIn,4);
%justafew=1:100:size(TrainDataIn,4);

Truth = TestDataOut(justafew,:);
%Truth = TrainDataOut(justafew,:);
est = zeros(size(Truth));

p=1;
for n = justafew
    TestInput = TestDataIn(:,:,:,n);
    %TestInput = TrainDataIn(:,:,:,n);
    tmp = predict(GWcalcNet , TestInput);
    % the output size will be (1 x Ncoils*Nnbrs*2)
    
    %{
    if (norm(tmp)==0)
        subplot(121)
        imagesc(TestInput)
        caxis([-1 1]*20)
        title('Bad')
    else
        subplot(122)
        imagesc(TestInput)
        caxis([-1 1]*20)
        title('good')
        drawnow
    end
    %}
    est(p,:) = tmp;
    p = p+1;
end
%

inds = find(est==0);
% frac = length(inds)/length(est(:))
%
% Truth(inds) = [];
% est(inds) =[];

figure

plot(Truth(:), est(:),'o')
rho = abs(corrcoef(Truth, est))
xlabel('Truth')
ylabel('Estimate')
legend(sprintf('Correlation %f', rho(2,1)))
line([min(abs(Truth)) max(abs(Truth))], [ min(abs(Truth)) max(abs(Truth))  ])




return

