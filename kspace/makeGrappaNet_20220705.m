function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
% function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
%
% first determine a ceter region for calibration assuming that this
% center region is well sampled.
%
% kx, ky, kz is the trajectory (cm*-1)
% data is the signal
% dim3 is num. pixels in each dimension of IMAGE space
% CalRadius is the fraction of k space used for calibration region
%
% This version learns how to calculate the weights by patches, instead of
% point-by-point
%

R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
Ncoils = size(data,2);

% the calibration region is determined by CalRadius relative to the Rmax
CalRadius = CalRadius * Kmax;  % convert from fraction to k space units.
CalRegion = find(R<CalRadius); % indices of data in the calibration region

dim3_cal = floor(dim3*CalRadius)+1;  % Now dim3 will be the num. pixels in cal region.
dim3_cal = 2*floor(dim3_cal/2) +1;   % Making sure it's an odd number

%%
fprintf('\rGRIDDING: Making cartesian calibration data from non cartesian data...\n')
% now interpolate spiral data into a cartesian grid
% will substitute this with kbgrid eventually
[subkx , subky, subkz ]= meshgrid( ...
    linspace(-CalRadius,CalRadius,dim3_cal(1)), ...
    linspace(-CalRadius,CalRadius,dim3_cal(2)), ...
    linspace(-CalRadius,CalRadius,dim3_cal(3)));


cal_grid = [];

for coilnum=1:Ncoils
    % calibration region for each coil
    % tmp = griddata(kx, ky, kz, data(:,coilnum), subkx, subky, subkz);
    fprintf('\rcoil  %d ... ', coilnum);
    tmp = grid3d_sinc(kx, ky, kz, data(:,coilnum), subkx, subky, subkz);
    
    cal_grid = [cal_grid tmp(:)];
end

% size of k-space voxels: distance between neighbors in the iso-tropic cartesian grid
dkx = abs(subky(2,1,1) - subky(1,1,1));

save workspace.mat
%%

fprintf('\nConstructing patches within calibration region ...\n ');

doRandomNeighbors = 1;

Nnbrs = 20;     % number of k space data locations for interpolation
% in each patch.  a patch is a constellation of points in
% kspace arround a location that we want to interpolate.

pklocs = [];    % The kspace coordinates of the data in each patch
% dimensions are: 3 x (Nequations-1)

Npatches = 1e4;  % number of times we'll calcuate the interpolation kernel


TrainDataIn =[];
TrainDataOut =[];

% vector with possible kspace targets (the middle half of the calibration
% region).  This is also the allowed kspace distance to possible neighbors
% in units of pixels (k-space)
%
center_coords = [-floor(dim3_cal(1)/16):floor(dim3_cal(1)/16)];
nbr_coords = [-floor(dim3_cal(1)*7/16):-1   1:floor(dim3_cal(1)*7/16)];
%nbr_coords = center_coords;
Nnbrs = 64;

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
% then we have to do it for all the missing data in all the coils

% make a vector to keep track of which coil is being used
coilvec =  [1:Ncoils]' / Ncoils;
coilvec = repmat(coilvec,1, Nnbrs)';
coilvec = coilvec(:);
coilvec = repmat(coilvec,Ncoils,1);

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

TrainDataIn = zeros(Nnbrs, 3,  Npatches);
TrainDataOut = zeros(Nnbrs*Ncoils*Ncoils, Npatches);

for p = 1:Npatches
    
    if doRandomNeighbors
        % option 1: the neigbors are randomized
        % k-neighborhood for interpolation
        % units of pixels (k-space)
        dx = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
        dy = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
        dz = nbr_coords(randi(length(nbr_coords), Nnbrs,1))';
    end
    
    
    Kneighbors(:) = 0;
    for n=1:Nequations
        % choose targets at random to serve as center of patch
        ind = randi(length(x), 1);
        %ind = n;
        
        patchcenter = [x(ind) y(ind) z(ind)];          % grid location of the target
        pklocs =  [x(ind)+dx  y(ind)+dy z(ind)+dz];    % grid locations of the neighbors
        
        % put all the coordintates into one matrix
        % the first k space location is the center of the patch
        pklocs = [patchcenter; pklocs];
        
        % the k space locations of interpolation neighbors are relative to
        % the center of the calibraiton region
        pklocs = pklocs + ceil(dim3_cal/2);
        
        % Extract the signals a those locations
        % (remember that the first one is the center/target)
        inds = sub2ind(dim3_cal, pklocs(:,1), pklocs(:,2), pklocs(:,3));
        
        tmp = cal_grid(inds(1),:);  % These are the target signals
        Ktargets(n,:) = tmp;
        
        tmp= cal_grid(inds(2:end),:);  % These are the neighbor signals (sources)
        Kneighbors(n,:) =  tmp(:)';
        
    end
    
    
    % Solve for the GRAPPA weights.
    [W, rnk, RMSE] = calc_grappa_kernel(Ktargets, Kneighbors);
    
    if mod(p,20)==0
        fprintf('\rSolving GRAPPA weights for patch ... %d \n', p);
        subplot(221)
        imagesc(abs(Kneighbors))
        title('neighbors')
        xlabel('Ncoil * Nneighbors')
        ylabel('Ntargets (Nequations)')
        
        
        subplot(222)
        imagesc(abs(Ktargets))
        title('targets')
        xlabel('Ncoil')
        ylabel('Ntargets (Nequations)')
        
        subplot(223)
        imagesc(abs(W))
        title('Weights')
        xlabel('Ncoil')
        ylabel('Ncoil x Nnbrs')
        
        subplot(224)
        plot(RMSE)
        ylabel('RMSE')
        xlabel('Ncoil')
        title(sprintf('rank: %d ', rnk))
        
        drawnow
    end
    
    W = W';
    
    % construct training data table
    % save the coordinates for this patch as a table
    % each entry (row) is of the form;
    % [ deltakx, deltaky deltakz coilnum ------> weight]
    %
    % in this looping order:
    %    for Npatches
    %       for Ncoils
    %           for Nneighbors
    %
    nbrLocs = [dx dy dz]*dkx;
    % neighbor locations RELATIVE to the target in this constelaltion/patch
    % units of cm^-1
    
    TrainDataIn(:,:,p) = nbrLocs;
    TrainDataOut(:,p) = W(:);
end

% shuffle anf spllit the training data into training and testing data sets.
Ntotal = size(TrainDataIn,3);
inds = randperm(Ntotal);

TestDataIn  = TrainDataIn(:,:, inds(1:Ntotal/10));
TrainDataIn = TrainDataIn(:,:, inds(Ntotal/10 +1:end));

TestDataOut  = TrainDataOut(:, inds(1:Ntotal/10))';
TrainDataOut = TrainDataOut(:, inds(Ntotal/10 +1:end))';

save TrainData.mat  TrainDataIn TrainDataOut TestDataOut TestDataIn
%%

% construct a neural net that will take in distances and coils and return the
% corresponding weights for an interpolation kernel
GWnet = [
    
imageInputLayer([Nnbrs 3]);


fullyConnectedLayer(Nnbrs*Ncoils );
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(Nnbrs*Ncoils*2);
batchNormalizationLayer;
reluLayer;

fullyConnectedLayer(Ncoils*Nnbrs*Ncoils);
regressionLayer
]

scale = norm((TrainDataOut(:)));
TrainDataOut =TrainDataOut(:)/scale;

% Arrange the training data into the right shape:
Nlearning = size(TrainDataIn,3);
Ntesting  = size(TestDataIn, 3);

TrainDataIn = reshape(TrainDataIn, Nnbrs, 3, 1, Nlearning);
TrainDataOut = abs(reshape(TrainDataOut, 1,1,Nnbrs*Ncoils*Ncoils, Nlearning));

TestDataIn = reshape(TestDataIn, Nnbrs, 3, 1, Nlearning);
TestDataOut = abs(reshape(TestDataOut, 1,1,Nnbrs*Ncoils*Ncoils, Nlearning));

fprintf('\nTraining Network ...\n ');
options = trainingOptions('adam',...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-4, ... %1e-4
    'WorkerLoad', ones([40 1]), ...
    'MaxEpochs',50, ...
    'LearnRateSchedule', 'piecewise', ...
    'Plots','training-progress');

[GWcalcNet, info] = trainNetwork(TrainDataIn, TrainDataOut , GWnet, options);

save GRAPPAnet.mat GWcalcNet info scale
%%

% scaling factor on output????

fprintf('\nTesting the network with data from calibration region...')
justafew=1:100;

%Truth = TestDataOut(justafew, :);
Truth = squeeze(TrainDataOut(:,:,justafew, :));
a = squeeze(TrainDataIn(:,:,:,justafew));

est = zeros(size(Truth));

for n = justafew
%    TestInput = TestDataIn(:,:,n);
    TestInput = a(:,:,n);
    tmp= predict(GWcalcNet , TestInput);
    est(n,:) = tmp * scale;
end

ERR = 100*mean(((est - Truth)./Truth).^2, 1);

figure
plot(abs(Truth(:)), abs(est(:)),'o')
rho = abs(corrcoef(Truth, est))
legend(sprintf('Correlation %f', rho(2,1)))
line([min(abs(Truth)) max(abs(Truth))], [ min(abs(Truth)) max(abs(Truth))  ])
%%

return

