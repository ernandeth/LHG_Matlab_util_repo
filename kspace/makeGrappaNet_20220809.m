function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
% function GWcalcNet = makeGrappaNet (kx, ky, kz, data, dim3, CalRadius)
%
% this version tried to convert neighbor locations into weights.
%
% THIS DOES NOT GRID the calibration data into a Cartesian space first
% it calculates grappa weights from the data in the central region
% without making it cartesian first.
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

R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R);
Ncoils = size(data,2);


CalRadius = CalRadius * Kmax;  % convert from fraction to k space units.
CalRegion = find(R<CalRadius ); % indices of data in the calibration region

% The coordinates within the calibration region 
ckx = kx(CalRegion);
cky = ky(CalRegion);
ckz = kz(CalRegion);
cdata = data(CalRegion, :);
cR = R(CalRegion);

% eliminate redundant locations
%{
Ndata = size(ckx,1);
repeats = [];
for n=2:Ndata
    if ckx(n) == ckx(n-1)
        if cky(n) == cky(n-1)
            if ckz(n) == ckz(n-1)
                
                repeats = [repeats;  n];
            end
        end
    end
end
cR(repeats) = [];
ckx(repeats) = [];
cky(repeats) = [];
ckz(repeats) = [];
cdata(repeats, :) = [];
%}
acs_coords = [ckx cky ckz];

save workspace.mat
%%

fprintf('\nConstructing patches within calibration region ...\n ');

doRandomNeighbors = 1;

Nnbrs = 15;     % number of k space data locations for interpolation
% in each patch.  a patch is a constellation of points in
% kspace arround a location that we want to interpolate.

nbr_coords = [];    % The kspace coordinates of the data in each patch
% dimensions are: 3 x (Nequations-1)

Npatches = 1000;  % number of times we'll calcuate the interpolation kernel

% vector with possible kspace targets (the middle half of the calibration
% region).  This is also the allowed kspace distance to possible neighbors
% in units of pixels (k-space)
%
target_inds = find(cR < CalRadius/4);

%  Targets: coordinates for the possible centers of the patches
tgtx = ckx(target_inds);
tgty = cky(target_inds);
tgtz = ckz(target_inds);


Nequations = round(Ncoils*Nnbrs*1.5);
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
% TrainDataOut = zeros( Nequations*Npatches, Ncoils*Nnbrs*2, Ncoils);
% we're going to do only one coil at a time
TrainDataOut = zeros( Nequations*Npatches, Ncoils*Nnbrs*2 );

train_cnt = 1;

for p = 1:Npatches
    
    
    Kneighbors(:) = 0;
    
    % COnstruct system of equations to solve for GRAPPA weights
    
    % First choose a set of RELATIVE distances for the neightbors
    nbr_coords = (rand(Nnbrs,3) * CalRadius - CalRadius/2)/2 ;
    
    for n=1:Nequations
        
        % choose target locations near the center of the calibration region
        % at random to serve as center of patch
        ind = randi(length(tgtx), 1);        
        tgt_coords = [tgtx(ind) tgty(ind) tgtz(ind)];          
        tgt = cdata(ind, :);  % These are the target complex signals
                                    % dimensions: 1 x Ncoils,
        Ktargets(n,:) = tgt;        % ... we stack them into Nequations x Ncoils
        
        %{
        %----------------------------
        % now choose some  random neighbors, 
        dist = (acs_coords - repmat(tgt_coords, size(acs_coords,1) ,1)).^2;
        dist = sqrt(sum(dist,2));
      
        % but not further than CalRadius/2 from the target
        inds = find( dist<CalRadius/2  & dist>0.1);
        inds = shuffle(inds);
        inds = inds(1:Nnbrs);
        nbr_coords = acs_coords(inds,:);
        % location of neighbors relative to the target
        % nbr_coords = nbr_coords - repmat(tgt_coords, Nnbrs,1);        
%----------------------------
        %}
        
        % for each neighbor location:
        % identify which k-space loc in the acquired trajectory is closest
        % to that specified neighbor location
        inds = zeros(Nnbrs,1);
        for b=1:Nnbrs
            % neightbor's positions are relative to the target
            nbrs = nbr_coords(b,:) + tgt_coords;
            
            nbrs = repmat(nbrs, size(acs_coords,1),1);
            dist = sqrt(sum((acs_coords - nbrs).^2 ,2 ));
            inds(b) = find(dist==min(dist));
        end
        
        % Neighbor complex signals (Nnbrs x Ncoils)
        tmp = cdata(inds,:);
        Kneighbors(n,:) =  tmp(:)';
        %{
        % plot the coordinates of the target and the neightbors
         
        plot3(ckx(inds), cky(inds), ckz(inds), 'ko')
        hold on
        plot3(tgt_coords(:,1), tgt_coords(:,2), tgt_coords(:,3),'r*')
        hold off
        axis([-1 1 -1 1 -1 1]*CalRadius)
        title(sprintf('Equation %d',n))
        drawnow
        %}

        nbrlocs = repmat(nbr_coords, Ncoils,1);

        % store relative locations for training the neural net
        TrainDataIn( :, :, 1, train_cnt) = nbrlocs;
        
        train_cnt = train_cnt +1;
        
    end
    
    
    % Solve for the GRAPPA weights.
    fprintf('\rSolving GRAPPA weights for patch ... %d ', p);
    %
    [W, rnk, RMSE] = calc_grappa_kernel(Ktargets, Kneighbors);
    
    % store the error for coil 1
    rmse(p) = RMSE(1);
    
    % Output training data will be the GRAPPA weights for each coil
    % we note that each equation in the block,
    % the answer was the same set of weights
    % We'll start trying to do just one coil at a time
    coilnum = 1;
    Wblock = repmat(W(:,coilnum)', Nequations,1);
    Wblock = [real(Wblock) imag(Wblock)];
    
    % weed out training data if the GRAPPA fit is not very good
    if mean(RMSE) > 5e4
        Wblock(:) = nan;
        fprintf('RMSE =%f too high - discard', mean(RMSE))
    end
    
    TrainDataOut(Nequations*(p-1)+1:Nequations*p,:) =  Wblock;
    %
    if 1 % mod(p,100)==0
        subplot(221)
        imagesc(abs(Kneighbors))
        colorbar
        title('neighbors')
        xlabel('N coils * N neighbors')
        ylabel('N equations')
        
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
        ylabel('N equations')
        
        
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

% fullyConnectedLayer(Ncoils);
% batchNormalizationLayer;
% reluLayer;
%
% fullyConnectedLayer(2*Nnbrs*Ncoils);
% batchNormalizationLayer;
% reluLayer;

% one coil only
fullyConnectedLayer(2*Nnbrs*Ncoils);
regressionLayer
]

scale = max((TrainDataOut(:)));
scale = 1;


% TrainDataIn(:,5,:,:) = TrainDataIn(:,5,:,:)/scale;
% TrainDataIn(:,4,:,:) = TrainDataIn(:,4,:,:)/scale;
% junk = randn(size(TrainDataIn));
% junk(:,4:5, 1, :) = TrainDataIn(:,4:5,1,:);
% TrainDataIn = junk;
% TrainDataIn(:,1:3,:,:) = TrainDataIn(:,1:3,:,:) * 100;

TrainDataOut = TrainDataOut/scale;
TestDataOut = TestDataOut/scale;

fprintf('\nTraining Network ...\n ');
options = trainingOptions('sgdm',...
    'Shuffle','every-epoch', ...
    'InitialLearnRate', 1e-4, ... %1e-4
    'WorkerLoad', ones([40 1]), ...
    'MaxEpochs',10, ...
    'LearnRateSchedule', 'piecewise', ...
    'Plots','training-progress');  % 'training-progress'

[GWcalcNet, info] = trainNetwork(TrainDataIn, TrainDataOut, GWnet, options);

save GRAPPAnet.mat GWcalcNet info scale
%

fprintf('\nTesting the network with data from calibration region...')
justafew=1:100:size(TestDataIn,4);
%justafew=1:100:size(TrainDataIn,4);

%Truth = TestDataOut(justafew,:);
Truth = TrainDataOut(justafew,:);
est = zeros(size(Truth));

p=1;
for n = justafew
    %TestInput = TestDataIn(:,:,:,n);
    TestInput = TrainDataIn(:,:,:,n);
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
%%

return

