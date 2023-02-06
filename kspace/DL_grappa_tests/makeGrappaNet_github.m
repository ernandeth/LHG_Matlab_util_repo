function [allCoilNets  network_corr] = makeGrappaNet_github (kx, ky, kz, data, dim3, CalRadius)
% function allCoilNets = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
%
% (C) Luis Hernandez-Garcia @ University of Michigan
% hernan@umich.edu
% 
% Learn a function that will map kspace neighbor locations
% into interpolation weights.
% 1 - It grids the calibration data into a Cartesian space 
% 2 - calculates grappa weights from those data
% 3 - trains a network to solve for those weights , given relative locations of neighbors
% 4 - returns  a network for each coil
%
% INPUTS:
% kx, ky, kz is the trajectory in k-space (cm*-1)
% data is the signal : (Npts*Nechoes) x Ncoils
% dim3 is num. pixels in each dimension of IMAGE space
% CalRadius is the fraction of k space used for calibration region
%
% OUTPUT:
% returns a cell array with Ncoils, .  Each cell is a network for each coi.
% 
% derived from version 20220729 - that version gave the weights so that we
% could interpolate the targets at a single coil at a time.  
% This new version will try to train a network that calculates weights for
% all the coils at the same time.

% STEP 1: Set up a calibration region
% Determine a center region for calibration, assuming that this
% center region is well sampled.
% the calibration region is determined by CalRadius relative to the Rmax
R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
Ncoils = size(data,2);
% size of k-space voxels: (1/FOV) (distance between neighbors in the iso-tropic cartesian grid)
dkx = 2*Kmax/dim3(1);    

% Find the in the calibration region (center)
CalRegion = find(R<CalRadius); 

CalRegion_dim = floor(dim3*CalRadius)+1;  % dimesnions of cartesian calibration data matrix   
% We will Upsample the region:  Double the resolution in the calibration region - 
%  --> more voxels in cal region
CalRegion_dim = 2*CalRegion_dim;
dkx = dkx/2;  

CalRadius = CalRadius * Kmax;    % convert from fraction to k space units.

%%
figure
fprintf('\rGRIDDING: Making cartesian calibration data from non cartesian data...\n')
% now interpolate spiral data into a cartesian grid
% will substitute this with kbgrid eventually
[subkx , subky, subkz ]= meshgrid( ...
    linspace(-CalRadius,CalRadius,CalRegion_dim(1)), ...
    linspace(-CalRadius,CalRadius,CalRegion_dim(2)), ...
    linspace(-CalRadius,CalRadius,CalRegion_dim(3)));

% dimensions: (calregion_Npix x Ncoils)
cal_grid = [];
dcf = [];
im = zeros(size(subkx));

for coilnum=1:Ncoils
     
    % calibration region for each coil
    fprintf('\rcoil  %d ... ', coilnum);
    
    %{
    tic
    % grid the calibration dat into a cartesian grid
    [tmp dcf kernel] = grid3d_lhg(...
        kx, ky, kz, ...
        data(:,coilnum), ...
        subkx, subky, subkz, ...
        2, dcf);
    toc
    %}
    
    % Alternatively ... use griddata
    tmp = griddata(kx, ky, kz, data(:,coilnum), subkx, subky, subkz,'linear');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    %
    
    % testing images for recon
    im = im + (fftshift(abs(ifftn(tmp)))).^2;
    %fk = fftshift(abs(fftn(kernel)));
    
    subplot(2,1,1)
    lightbox(log(abs(tmp ))); 
    title(['coil ' num2str(coilnum) ' calibration region interpolated']); 
    
    subplot(2,1,2)
    lightbox(sqrt(im)); 
    title(['coil: ' num2str(coilnum) ' cal region in image domain']); 
    
%     subplot(3,1,3)
%     %lightbox(abs(fk(10:end-10,10:end-10, 10:end-10)),[],5,[]); 
%     title('kernel rolloff (image domain)');
    
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


save workspace.mat
%%

% STEP 2: GRAPPA make patches to interpolate within the calibration region.

% a known TARGET
fprintf('\nConstructing patches within calibration region ...\n ');
load workspace.mat

doRandomNeighbors = 1;

Nnbrs = 21;     % number of k space data locations for interpolation
Nnbrs = 35;
% in each patch.  a patch is a constellation of points in
% kspace arround a location that we want to interpolate.

nbrlocs = [];    % The kspace coordinates of the data in each patch
% dimensions are: 3 x (Nequations-1)

Npatches = 2000;  % number of times we'll calcuate the interpolation kernel
%Npatches = 1e4;

% Number of equations in a patch.
% The number of unknowns in the interpolation kernel
% is Ncoils*Nnbrs,
% so we should need AT LEAST Ncoils*Nbrs  equations
% to solve that many unknowns per grappa patch
Nequations = round(Ncoils*Nnbrs);
Nequations = 100;  % somehow this works better??


% vector with possible kspace targets (the middle half of the calibration
% region).  
center_coords = [-floor(CalRegion_dim(1)/4):floor(CalRegion_dim(1)/4)]

%  Targets: grid of possible coordinates for the center of the patches
[x y z] = meshgrid(center_coords, center_coords, center_coords);
x = x(:);
y = y(:);
z = z(:);

% Neighbors:  coordinates for the neighbors.  They are equally spaced
% on the neighborhood around the center
% (use the neighbors within the CalRegion_dim/8 region of the target)
% allowed kspace distance (voxel offset) to the possible neighbors
% in units of pixels (k-space)
% nbr_coords = [-floor(CalRegion_dim(1)/8):-1  1:floor(CalRegion_dim(1)/8)];
nbr_coords = [-8:-1 1:8]; % recall that the resolution has been doubled.

% Neighbors: grid of possible voxel offsets for the neighbors
[dx dy dz] = meshgrid(nbr_coords, nbr_coords, nbr_coords);
dx = dx(:);
dy = dy(:);
dz = dz(:);


% STEP 3: calculate the interpolation weights for the patches 
%
% GRAPPA interpolation does this:
%
% targetSignals = Weights .* nbrSignals
%
% nbrSignals    % dimensions: Nequations x Ncoils*Nnbrs
% Weights       % dimensions: Ncoils*Nnbrs x Ncoils
% targetSignals      % dimensions: Nequations x Ncoils
%

nbrSignals = zeros(Nequations, Ncoils*Nnbrs);
targetSignals = zeros(Nequations, Ncoils);

% Allocate space for the targets, neighbors, and weights that we will use
% for training and testing the network
% TrainDataIn = zeros( Ncoils*Nnbrs , 3, 1, Npatches*Nequations);
% BUT .... in the line above, the  coordinates repeat for every coil:
% --> Redundant information!  use only Nnbrs x 3
TrainDataIn = zeros( Nnbrs , 3, 1, Npatches*Nequations);

% The output of the network are  the weights of the kernel for each neighborhood
TrainDataOut = zeros(  Ncoils*Nnbrs*2 , Ncoils, Nequations*Npatches);
TS = size(TrainDataOut, 3);
badSamps = zeros(length(TS),1);

% Now start generating training data by solving the interpolation weights
% (GRAPPA kernel) for all the pathches 
train_cnt = 1;
figure
for p = 1:Npatches
    
    % neighbor locations: 
    % the neigbors are randomized for this patch
    % These are the offsets of the neighbor locations relative to the target
    % (still in units of voxels)
    dx = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    dy = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    dz = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    nbrlocs =  [dx dy dz];

   
    nbrSignals(:) = 0;
    
    % COnstruct system of equations to solve for GRAPPA weights
    for n=1:Nequations
        % choose one target at random to serve as center of patch (target)
        ind = randi(length(x), 3, 1);
        patchcenter = [x(ind(1)) y(ind(2)) z(ind(3))];  % grid location of the target
        nbrlocs =     patchcenter + [dx  dy dz];        % grid locations of the neighbors
       
        %{
        % ------Plots to make sure that we're in the right place
        cog = mean(nbrlocs,1);
        plot3(x(ind), y(ind), z(ind),'ro')
        hold on
        plot3(x(ind)+dx,  y(ind)+dy, z(ind)+dz,'kx')
        plot3(cog(1), cog(2), cog(3), 'bo')
        hold off
        axis([-1 1 -1 1 -1 1]*10)
        pause(0.1)
        % ----------------------------------------------------------
        %}
        
        % Extract the signals at those locations
        % targets:
        patchcenter = patchcenter + CalRegion_dim/2 +1;
        target_ind = sub2ind(CalRegion_dim, patchcenter(1), patchcenter(2), patchcenter(3));
        
        % Target complex signals (1 x Ncoils)
        tgt = cal_grid(target_ind,:);  % These are the target complex signals at each coil
        targetSignals(n,:) = tgt;   % ... we stack them into a matrix 
                                    % dims:  Nequations x Ncoils
        
        % neighbors: 
        nbrlocs = nbrlocs + CalRegion_dim/2 +1;
        nbr_inds = sub2ind(CalRegion_dim, nbrlocs(:,1), nbrlocs(:,2), nbrlocs(:,3));
        
        % Neighbor complex signals (1 x Nnbrs*Ncoils)
        nbrsig = cal_grid(nbr_inds,:);
        nbrSignals(n,:) =  nbrsig(:)';   % flattened into a single row
                                         % order: neighbors go together, then
                                         % the next coil ... 
                                         % ... we stack them into a matrix
                                         % dims:  Nequations x Ncoils*Nnbrs
                                         
        % We will use nbrlocs to make a table of locations
        % for training the input of the neural net later
        % dimension: (Nnbrs x 3 x 1 x NtrainingSamps)
        
        % IMPORTANT: positions are relative to the target !
        % Use kspace units , ie : (cm^-1), 
        rel_nbrlocs = (nbrlocs - patchcenter) *dkx ;
      

        TrainDataIn( :, :, 1, train_cnt) = rel_nbrlocs;
        train_cnt = train_cnt +1;
        
    end
    
    % Now solve the GRAPPA inverse for each patch.  We will store the
    % weights for training the output of the neural net
    %
    % Theory:  interpolate targets from neightbors
    %
    %    targetSignals = nbrSignals * Weights
    %
    % dimensions:
    % targetSignals :    Nequations x Ncoils
    % nbrSignals :  Nequations x (Nnbrs*Ncoils)
    % W:            (Nnbrs*Ncoils) x Ncoils
    %     
    % Invert this equation to solve for the weights
    % (the weigths are the same for each block of equations) 
    % 
    fprintf('\rSolving GRAPPA weights for patch ... %d ', p);
    
    [W, rnk, RMSE] = calc_grappa_kernel(targetSignals, nbrSignals);

    % store the error for coil 1
    % ... but weed out training data if the GRAPPA fit is not good
    rmse(p) = RMSE(1);
    fprintf('RMSE =%f ', mean(RMSE));
    if mean(RMSE) > 1e-8
        W(:) = nan;
        fprintf('RMSE too high - discard')
        badSamps(Nequations*(p-1)+1:Nequations*p) = 1;
    end
    
    %{
    hold off; plot(abs(W(:,1))/max(abs(W(:,1))),'r')
    hold on; plot(abs(nbrSignals(1,:))/max(abs(nbrSignals(1,:))),'b')
    drawnow
    %}
    
    % note the resulting W are complex weights: so we will store as real,
    % imag for the output of the neural net.
    % Later. The output of the Neural net will be 
    % the GRAPPA weights for each neighbor and each coil 
    % needed to interpolate a single target point in a single coil
    % from all the neightbors in all the coils.
    % Dimension: (Nnbrs*Ncoils*2) x 1 ...
    W2 = [real(W); imag(W)];

  
    % W is the solution to the whole system of equations that we constructed,
    % so each of these equations can be used as a training example
    %  -->    (Nnbrs*Ncoils*2) x Ncoils x Nequations
    Wblock= repmat(W2, 1,1,Nequations);
    
    % The out put of the network is the weights. These weights are the
    % output corresponding to each ot the above inputs
    TrainDataOut(:,:, Nequations*(p-1)+1:Nequations*p) =  Wblock;
    
    if mod(p,200)==0
        subplot(221)
        imagesc(abs(nbrSignals))
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
        imagesc(abs(targetSignals))
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
figure

%clear nbrSignals targetSignals Wblock

% remove Training data containing NaNs
%
inds = find(badSamps);
TrainDataIn(:,:,:,inds) = [];
TrainDataOut(:,:,inds) = [];
%}

% split the training data into training and testing data sets.
Ntrain = size(TrainDataIn,4);
%
% Grab a fraction of them for testing
inds = randperm(Ntrain, ceil(Ntrain/20));
TestDataIn = TrainDataIn( :,:,1 , inds);
TestDataOut = TrainDataOut(:,:,inds);

% ... and remove those from the training data set
TrainDataOut(:,:, inds) = [];
TrainDataIn( :,:,:, inds) = [];

fprintf('\nSaving Training Data ...');
whos Test* Train*
%save -v7.3 TrainData.mat Test* Train* N*
%%
% Goal: Learn a function that will map relative neighbor locations into complex
% weights for interpolation
%
% how: construct a neural net that will map the relative kspace locations into
% the relative interpolation weights for one coil at a time.
% input dimensions:  (Nnbrs x 3)
% output dimensions:  (Nnbrs*Ncoils*2 x 1)
%%
%load TrainData.mat
%
GWnet = [
    
    imageInputLayer([Nnbrs 3]);
    
    fullyConnectedLayer(Ncoils*Nnbrs*2);
    batchNormalizationLayer;
    reluLayer;
    
    fullyConnectedLayer(Ncoils*Nnbrs*10);
    batchNormalizationLayer;
    reluLayer;

    fullyConnectedLayer(Ncoils*Nnbrs*10);
    batchNormalizationLayer;
    reluLayer;

    
    fullyConnectedLayer([Nnbrs*Ncoils*2]);
    regressionLayer
]
%}
scale = max((TrainDataOut(:)));
scale = 1;

%TrainDataOut = TrainDataOut;
%TestDataOut = TestDataOut;

fprintf('\nTraining Networks ...\n ');
options = trainingOptions('sgdm',...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-3, ... %1e-4
    'WorkerLoad', ones([40 1]), ...
    'MaxEpochs',10, ...
    'LearnRateSchedule', 'piecewise', ...
    'Plots','none'); % 'training-progress');

% Now train this network once for each coil
for n=1:Ncoils
   TDO = squeeze(TrainDataOut(:,n,:))';
        
   [coilNet, info] = trainNetwork(...
        TrainDataIn, ...
        TDO, ...
        GWnet, ...
        options);
   
    allCoilNets{n} = coilNet;
end

%fprintf('\nSaving networks ...')
%save GRAPPAnet.mat  allCoilNets


%
figure
fprintf('\nTesting the network with data from calibration region...')
network_corr = zeros(Ncoils,1);
for coilnum = 1:Ncoils
    justafew=1:100:size(TestDataIn,4);
            
    TestInput = TestDataIn(:,:,:,justafew);
    Truth = squeeze(TestDataOut(:,coilnum,justafew));
    est = zeros(size(Truth));
    
    p=1;
    for n =1:length(justafew)
        tmp = predict(allCoilNets{coilnum} , TestInput(:,:,:,n));
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
        est(:,p) = tmp;
        p = p+1;
    end
    %
    
    inds = find(est==0);

    
    %
    figure
    plot(Truth(:), est(:),'o')
    rho = abs(corrcoef(Truth, est))
    xlabel('Truth')
    ylabel('Estimate')
    legend(sprintf('Correlation %f', rho(2,1)))
    line([min(abs(Truth)) max(abs(Truth))], [ min(abs(Truth)) max(abs(Truth))  ])
    title(['Coil ' num2str(coilnum)])
    pause(0.1)
    
    %}
    network_corr(coilnum) = rho(2,1);
end

return
