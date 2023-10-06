function [allCoilNets  network_corr] = makeGrappaNet (kx, ky, kz, data, dim3, CalFraction)
% function allCoilNets = makeGrappaNet (kx, ky, kz, data, dim3, CalFraction)
%
% (C) Luis Hernandez-Garcia @ University of Michigan
% hernan@umich.edu
% 
% Learn a function that will map kspace  locations
% into k space data 
%
% INPUTS:
% kx, ky, kz is the trajectory in k-space (cm*-1)
% data is the signal : (Npoints_echo*Nechoes) x Ncoils
% dim3 is num. pixels in each dimension of IMAGE space
% CalFraction is the fraction of k space used for calibration region
%
% OUTPUT:
% returns a cell array.  Each cell is a network for each coi.
% 
%

% 
% Skipping this section and reading whatever was saved in the
% workspace.mat :

% figure out how many coisl are present
Ncoils = size(data,2);

% STEP 1: Set up a calibration region
% Determine a center region for calibration, assuming that this
% center region is well sampled.
% the calibration region is determined by CalFraction relative to the Rmax
R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
CalRegion = find(R < CalFraction*Kmax);   % Find the points that belong in the calibration region (center)

dkx = Kmax/dim3(1);    % size of k-space voxels: 
                         %  (1/FOV) (distance between neighbors in the iso-tropic cartesian grid)

%  Upsample:  Double the resolution in the calibration region - 
%  --> more voxels in cal region
% CalRegion_dim = 2*CalRegion_dim;
% dkx = dkx/2;  


fprintf('\rGRIDDING: Making cartesian calibration data from non cartesian data...\n')
% now interpolate spiral data into a cartesian grid.  The GRAPPA
% interpolation will be done on these data.
% we have to do this for each coil
cal_grid = [];
dcf = [];
sampFactor = 1
% create the mesh for interpolation
CalKmax = CalFraction * Kmax;    % convert from fraction to k space units.
CalRegion_dim = sampFactor* floor(dim3*CalFraction);  % size of cartesian calibration data matrix   
dkx = dkx/sampFactor;  

[subkx , subky, subkz ]= meshgrid( ...
    linspace(-CalKmax,CalKmax,CalRegion_dim(1)), ...
    linspace(-CalKmax,CalKmax,CalRegion_dim(2)), ...
    linspace(-CalKmax,CalKmax,CalRegion_dim(3)));
im = zeros(size(subkx));

% Make a Gaussian filter in k-space.  WIll use this for 
% smoothness as regularization.
% Gfilter = exp(-(subkx.^2 + subky.^2 + subkz.^2) ./ CalRegion_dim(1).^2 );

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
    
    % apply a filter to the data -  smoothness in image domain
    % this should act as a regularizer of sorts.
    % tmp = tmp .* Gfilter;
    % tmp = smooth3(tmp, 'gaussian',5);  <---- this never seems to help

    % testing images for recon
    im = im + (fftshift(abs(ifftn(tmp)))).^2;
    %fk = fftshift(abs(fftn(kernel)));
    
    subplot(2,1,1)
    lightbox(log(abs(tmp ))); 
    title(['coil ' num2str(coilnum) ' calibration region (gridded)']); 
    
    subplot(2,1,2)
    lightbox(sqrt(im)); 
    title(['coil: ' num2str(coilnum) ' cal region in image domain']); 
    
%     subplot(3,1,3)
%     %lightbox(abs(fk(10:end-10,10:end-10, 10:end-10)),[],5,[]); 
%     title('kernel rolloff (image domain)');
    
    drawnow
    %
    cal_grid = [cal_grid tmp(:)];
end
clear tmp
% create a calibration region from the original data
% it will consist of the center fraction (CalFraction) of the sampled data
% and it will be on a cartesian grid.
%
% We will sampFactor the k-space calibration region - 
%  --> more voxels in cal region
%{
sampFactor = 2;
CalRegion_dim = sampFactor*CalRegion_dim;
dkx = dkx/sampFactor;  

cal_grid = grid2cartesian(kx, ky, kz, data, dim3*sampFactor, CalFraction);
%}

%save workspace.mat

% in case we want to skip this ....
%}

%} end the gridding sction

%%
% CREATE TRAINING DATA
% STEP 2: GRAPPA make patches to interpolate within the calibration region.

% a known TARGET
fprintf('\nConstructing patches within calibration region ...\n ');
%load workspace.mat

doRandomNeighbors = 1;
     
% number of k space data neighbors for interpolation
% in each patch.  a patch is a constellation of points in
% kspace arround a location that we want to interpolate.
%Nnbrs = 35;  <--------LHG 11.17.22 - this worked for the digital phantom
Nnbrs = 7;
Nnbrs = 10;
Nnbrs = 15;
%Nnbrs = 35;


nbrlocs = [];    % The kspace coordinates of the data in each patch
% dimensions are: 3 x (Nequations-1)

Npatches = 1000;  % number of times we'll calcuate the interpolation kernel
Npatches = 2000;  % number of times we'll calcuate the interpolation kernel
Npatches = 3000; % 10/01/23) Next time do 5000
Npatches = 5000; % 10/02/23 improved the slopes, but not the correlations?

% Number of equations in a patch.
% The number of unknowns in the interpolation kernel
% is Ncoils*Nnbrs,
% so we should need AT LEAST Ncoils*Nbrs  equations
% to solve that many unknowns per grappa patch
Nequations = round(Ncoils*Nnbrs);
Nequations = 100;  % somehow this works better??
Nequations = round(Ncoils*Nnbrs)/2; % this is 160 for 10 neighbors
Nequations = round(Ncoils*Nnbrs);


% vector with possible kspace targets (the middle half of the calibration
% region).  
center_coords = [-floor(CalRegion_dim(1)/4):floor(CalRegion_dim(1)/4)];

%  Targets: grid of possible coordinates for the center of the patches
[x y z] = meshgrid(center_coords, center_coords, center_coords);
x = x(:);
y = y(:);
z = z(:);

% Neighbors: relative coordinates for the neighbors.  They are on a cartesian grid
% in the neighborhood of the target - in the center calibration region.
% Allowed kspace distance (voxel offset) to the possible neighbors
% in units of pixels (k-space)
nbr_coords = [-6:6]; % recall that the resolution has been doubled.


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

for p = 1:Npatches
    
    % neighbor locations: 
    % the neigbors are randomized for this patch
    % These are the offsets of the neighbor locations relative to the target
    % (still in units of voxels)
    dx = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    dy = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    dz = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    nbrlocs =  [dx dy dz];

    % reject patches where one of the points is at 0,0,0
    for n=1:size(nbrlocs,1)
        if nbrlocs(n,:)==[0,0,0]
            nbrlocs(n,:)=shuffle([1,0,0]);
        end
    end

   
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
        plot3(patchcenter(1), patchcenter(2), patchcenter(3),'ro')
        hold on
        plot3(nbrlocs(:,1), nbrlocs(:,2) , nbrlocs(:,3),'kx')
        plot3(cog(1), cog(2), cog(3), 'bo')
        hold off
        axis([-1 1 -1 1 -1 1]*dim3(1)*CalFraction)
        drawnow
        %pause(0.1)
        % ----------------------------------------------------------
        %}
        
        % Extract the signals at those locations
        % targets:
        patchcenter = round(patchcenter + CalRegion_dim/2 );
        target_ind = sub2ind(CalRegion_dim, patchcenter(1), patchcenter(2), patchcenter(3));
        
        % Target complex signals (1 x Ncoils)
        tgt = cal_grid(target_ind,:);  % These are the target complex signals at each coil
        targetSignals(n,:) = tgt;   % ... we stack them into a matrix 
                                    % dims:  Nequations x Ncoils
        
        % neighbors' coordinates: 
        nbrlocs = round(nbrlocs + CalRegion_dim/2) +1 ;
        nbr_inds = sub2ind(CalRegion_dim, nbrlocs(:,1), nbrlocs(:,2), nbrlocs(:,3));
        
        % Neighbor's complex signals (1 x Nnbrs*Ncoils)
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

    % fprintf('\rSolving GRAPPA weights for patch ... %d ', p);
    
    [W, rnk, RMSE] = calc_grappa_kernel(targetSignals, nbrSignals);

    % store the error for coil 1
    % ... but weed out training data if the GRAPPA fit is not good
    rmse(p) = RMSE(1);
    fprintf('% RMSE =%f ', mean(RMSE));
    if mean(RMSE) > 1e-8
        W(:) = nan;
        fprintf('\rSolving GRAPPA weights for patch ... %d ->', p);
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
        figure (21)
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
        title(sprintf('Error(rnk=%d)', rnk))
        xlabel('coils')
        ylabel('NRMSE (%)')
        
        drawnow
    end
    %}
    
end


%clear nbrSignals targetSignals Wblock

% remove Training data containing NaNs
%
fprintf('\nRemoving BAD training samples ...\n ');
inds = find(badSamps);
TrainDataIn(:,:,:,inds) = [];
TrainDataOut(:,:,inds) = [];
%}

%{ 
% randomize the order of the training data
fprintf('\nRandomizing order of training samples ...\n ');
shuffle_inds = randperm(size(TrainDataIn,4)); 
TrainDataIn = TrainDataIn(:,:,:,shuffle_inds);
TrainDataOut = TrainDataOut(:,:,shuffle_inds);
%}
% split the training data into training and testing data sets.
fprintf('\nSplitting training data into training and testing ...\n ');
Ntrain = size(TrainDataIn,4);
%
% Grab a fraction of them for testing
inds = randperm(Ntrain, ceil(Ntrain/20));
TestDataIn = TrainDataIn( :,:,1 , inds);
TestDataOut = TrainDataOut(:,:,inds);

% ... and remove those from the training data set
notinds = ones(Ntrain,1);
notinds(inds) = 0;
notinds = find(notinds);

%{
fprintf('\nSaving Training Data ...');
whos Test* Train*
save -v7.3 TrainData.mat Test* Train* N*
fprintf('done');
%}
%
% CREATE NETWORK 
%
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
% - good one 10/01/23
GWnet = [
    
    imageInputLayer([Nnbrs 3]);
    
    fullyConnectedLayer(Nnbrs*3);
    batchNormalizationLayer;
    reluLayer;


    fullyConnectedLayer(Ncoils*Nnbrs*25);
    batchNormalizationLayer;
    reluLayer;

    fullyConnectedLayer(Ncoils*Nnbrs*25);
    batchNormalizationLayer;
    reluLayer;


    fullyConnectedLayer([Nnbrs*Ncoils*2]);
    regressionLayer
]
%}

%-------------------------------------------
% experimental :great correlations during testing!!!! 10.02.2023
% but some coils have slope >> 1
GWnet = [
    
    imageInputLayer([Nnbrs 3]);
    
    fullyConnectedLayer(Nnbrs*3);
    batchNormalizationLayer;
    reluLayer;


    fullyConnectedLayer(Ncoils*Nnbrs*50);
    batchNormalizationLayer;
    reluLayer;

    fullyConnectedLayer([Nnbrs*Ncoils*2]);
    regressionLayer
]
%-------------------------------------------
% Now trying to make it wider ... 10/02/2023
% not very good with 3000 patches on some coils.  Great in others.
% trying now with 5000 patches - starting 10/02/2023. will check tomorrow ... 
GWnet = [
    
    imageInputLayer([Nnbrs 3]);
    
    fullyConnectedLayer(Nnbrs*3);
    batchNormalizationLayer;
    reluLayer;


    fullyConnectedLayer(Ncoils*Nnbrs*100);
    batchNormalizationLayer;
    reluLayer;

    fullyConnectedLayer([Nnbrs*Ncoils*2]);
    regressionLayer
]
%-------------------------------------------
scale = max((TrainDataOut(:)));
scale = 1;

%TrainDataOut = TrainDataOut;
%TestDataOut = TestDataOut;

% TRAIN NETWORK
%
fprintf('\nTraining Networks ...\n ');
options = trainingOptions('sgdm',...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-4, ... %1e-4
    'WorkerLoad', ones([40 1]), ...
    'MiniBatchSize', 20, ...  % reduce this (20) to avoid running out of memory (11.18.22)
    'MaxEpochs',10, ... % used to be 10
    'LearnRateSchedule', 'piecewise', ...
    'Plots', 'training-progress');
    
network_corr = zeros(Ncoils,1);

% Now train this network once for each coil
for coilnum=1 :Ncoils
   TDO = squeeze(TrainDataOut(:,coilnum,:))';
        
   [coilNet, info] = trainNetwork(...
        TrainDataIn, ...
        TDO, ...
        GWnet, ...
        options);
   
    allCoilNets{coilnum} = coilNet;
    
    

    %  TEST NETWORK for this coil
    fprintf('\nTesting the network with data from calibration region...')
    gpuDevice(1)
    justafew=1:100:size(TestDataIn,4);
            
    %TestInput = TestDataIn(:,:,:,justafew);
    %Truth = squeeze(TestDataOut(:,coilnum,justafew));

    % a sanity check
    TestInput = TrainDataIn(:,:,:,justafew);
    Truth = squeeze(TrainDataOut(:,coilnum,justafew));

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

    rho = corrcoef(Truth(:), est(:));
    network_corr(coilnum) = rho(2,1);

    figure
    plot(Truth(:), est(:),'.')
    xlabel('Truth')
    ylabel('Estimate')
    legend(sprintf('Correlation %f', rho(2,1)))
    line([min(abs(Truth(:))) max(abs(Truth(:)))], [ min(abs(Truth(:))) max(abs(Truth(:)))  ])
    title(['Coil ' num2str(coilnum)])
    drawnow
    %}
end
%%
%fprintf('\nSaving networks ...')
    %save GRAPPAnet.mat  allCoilNets

return

