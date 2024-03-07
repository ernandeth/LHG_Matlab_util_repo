function [netParms rho] = makeGrappaNet (kx, ky, kz, data, dim3, CalFraction)
% function netParms = makeGrappaNet (kx, ky, kz, data, dim3, CalFraction)
%
% (C) Luis Hernandez-Garcia @ University of Michigan
% hernan@umich.edu
%
% This function will train a neural network to
% learn the GRAPPA interpolation kernel weights
% from arbitrary neighbor locations relative to the target
%
% 2023.01.02 The training is done in non-cartesian data
% 2023.12.02 The training can  now be done in Cartesian data
%
% INPUTS:
% kx, ky, kz is the trajectory in k-space (cm*-1)
% data is the signal : (Npoints_echo*Nechoes) x Ncoils
% dim3 is num. pixels in each dimension of IMAGE space
% CalFraction is the fraction of k space used for calibration region
%
% OUTPUT:
% netParms: a structure containing the network 
% rho: correlation between testing data and truth

%
% Skipping this section and reading whatever was saved in the
% workspace.mat :

train_cartesian=0

% figure out how many coils are present
Ncoils = size(data,2);

% STEP 1: Set up a calibration region
% Determine a center region for calibration, assuming that this
% center region is well sampled.
% the calibration region is determined by CalFraction relative to the Rcal
R = sqrt(kx.^2 + ky.^2 + kz.^2);
Kmax = max(R)
dkx = 2*Kmax/dim3(1)    % size of k-space voxels:

Rinterp = 3*dkx                  % Radius of interpolation neighborhood
Kmax_src = CalFraction * Kmax + Rinterp    % Radius of calibration region (sources)
Kmax_tgt = CalFraction * Kmax ;     % Radius of calibration regions(targets)
if train_cartesian
    Kmax_tgt = Kmax_tgt/sqrt(2);
end

% remove duplicates of the center of k-space
kcenter_inds= find(R < 1e-5);
%first average all the centers of k space for each coil
for c=1:Ncoils   
    data(kcenter_inds, c) = mean(data(kcenter_inds,c)) ;
end
% then remove all the center points except for the first one in each coil
kcenter_inds = kcenter_inds(2:end);
data(kcenter_inds,:)=[];
% and remove the coordinates too (0,0,0)
kx(kcenter_inds) = [];
ky(kcenter_inds) = [];
kz(kcenter_inds) = [];
R = sqrt(kx.^2 + ky.^2 + kz.^2);
% adjust sizes    
Nsamples = size(data,1);

fprintf('\nTry griddata first ...\n')
% Interpolate from non-cartesian to cartesian grid using griddata
% do it coil by coil - note that each coil is done indepently of the
% others.  These pixel locations will be the targets for the Network
% interpolation

% dimensions of the cartesian grid for the calibration region
dim_cal = round(dim3 *CalFraction);
Npix_cart = dim_cal(1) * dim_cal(2) *dim_cal(3);
interped_data_gd = zeros(Npix_cart , Ncoils);

% Make the Cartesian grid
[xc , yc, zc ]= meshgrid( ...
    linspace(-Kmax_tgt,Kmax_tgt,dim_cal(1)), ...
    linspace(-Kmax_tgt,Kmax_tgt,dim_cal(2)), ...
    linspace(-Kmax_tgt,Kmax_tgt,dim_cal(3)));


parfor coilnum=1:Ncoils
    fprintf('\nInterpolating Cartesian Grid from coil %d using griddata',  coilnum)
    tmp = griddata(kx, ky,kz , data(:,coilnum), xc, yc, zc,'natural');
    tmp(isnan(tmp)) = eps;
    tmp(isinf(tmp)) = eps;
    interped_data_gd(:,coilnum) = tmp(:);
end

%%
% CREATE TRAINING DATA SET
fprintf('\nConstructing patches within calibration region ...\n ');
%load workspace.mat
doRandomNeighbors = 1;

% Find the possible target points that belong in the calibration region (center)
% the targets are a random collection Cartesian points
% in the middle half of the calibration region (Cal_tgt_inds)
if train_cartesian
    Cal_tgt_inds=1:Npix_cart;
else
    Cal_tgt_inds = find(R < Kmax_tgt );   
end

% Find the possible source points that belong in the calibration region (center)
% The neighbor sources will be non-Cartesian points within Rinterp of the target
% So the number of samples is reduced to the calibration region only
Cal_src_inds = find(R < Kmax_src );   
Nsamples = length(Cal_src_inds);

% make some convenient vectors for the k-space coordinates
kc = [xc(:) yc(:) zc(:)];  % cartesian grid coordinates
ks = [kx ky kz];           % experimental coordinates (non-cartesian) 

% reduce the data set to the training region
ks = ks(Cal_src_inds,:);
cal_data = data(Cal_src_inds,:);

% A "patch" is a constellation of points in
% kspace around a target location that we want to interpolate.
Npatches = 25000;  % number of times we'll calculate the interpolation kernel

% Number of k space source neighbors for interpolation
% in each patch. 
Nnbrs = 5;


% allocate space
dist = zeros(Nsamples,1);
dens = zeros(1,Npatches);
dist = zeros(1,Nsamples);
nbr_locs = zeros(Npatches, Nnbrs,3); % The kspace coordinates of the sources in each patch
Sdat = zeros(Npatches, Nnbrs*Ncoils);
Tdat = zeros(Npatches, Ncoils);


% Make a list of the targets at random order
% we need Npatches from the Cal_tgt_inds 
% can repeat them
idx = randi(length(Cal_tgt_inds),1, Npatches);
Cal_tgt_inds =  idx;

badPatch=zeros(1,Npatches);

% Loop for a number of target points in k-space data
% p is the index the target location to interpolate from neighbors
% pn is the counter of each patch
for pn=1:Npatches

    p = Cal_tgt_inds(pn);

    % Calc distance from p to all the other points in the
    % trajectory inside the Calibration Region (source)
    
    if train_cartesian
        % (#Cartesian Training)
        dist = sqrt(sum( (kc(p,:)-ks).^2 ,2));  % cartesian targets
    else
        dist = sqrt(sum( (ks(p,:)- ks).^2 ,2)); % spiral targets
    end
    
    % Find which points are within the Rinterp distance
    % Rinterp is the neighborhood for interpolation
    nbr_inds = find(dist < Rinterp  &  dist > 0.01*dkx);   
    dens(pn) = length(nbr_inds);      % save the sampling density for later use

    if(length(nbr_inds)<Nnbrs)
        badPatch(pn)=1;
        %{
        title(sprintf('pix %d Only %d neighbors !',p, length(nbr_inds)))

        plot3(ks(nbr_inds,1) , ks(nbr_inds,2), ks(nbr_inds,3),'r.'); hold on
        % plot3(ks(p,1), ks(p,2), ks(p,3), 'b.'); hold on
        % (#Cartesian Training)
        plot3(xc(p), yc(p), zc(p), 'b.');
        axis([-1 1 -1 1 -1 1]*Kmax)
        drawnow
        %}
    else
        % Choose a random set of Nnbrs in that Kmax_src region
        if doRandomNeighbors==1
            tmp = randperm(length(nbr_inds), Nnbrs);
            inds = nbr_inds(tmp);
        else
            % ... OR ... Choose the closest points
            [dist inds] = sort(dist);
            inds = inds(1:Nnbrs);
        end

        if mod(pn, 500)==0

            %figure(1)
            fprintf('\nMaking training data from calibration region patch num. %d \n ', pn);

            plot3(ks(inds,1) , ks(inds,2), ks(inds,3),'g.'); hold on
            if train_cartesian
                plot3(xc(p), yc(p), zc(p), 'r.'); hold on
            else
                plot3(ks(p,1), ks(p,2), ks(p,3), 'r.'); hold on
            end

            axis([-1 1 -1 1 -1 1]*Kmax_src)

            drawnow
            % figure(2)
            % tmp=nbr_locs(1:pn,:,:);
            % hist(tmp(:), 1000)
            % drawnow
        end

        % Now calculate locations of neighbors RELATIVE to the target
        % for this patch - this will be the INPUT to the network
        % the output of the network should be the target data 

        % Neighbor locations:   size (Nnbrs x 3)
        if train_cartesian
            nbr_locs(pn,:,:) = [ ks(inds,1)-xc(p) ,  ks(inds,2)-yc(p) ,   ks(inds,3)-zc(p)];
        else
            nbr_locs(pn,:,:) = [ ks(inds,1)-ks(p,1) ,  ks(inds,2)-ks(p,2) ,   ks(inds,3)-ks(p,3)];
        end
      

        % The output of the network is the target signal (Tdat):
        if train_cartesian
            Tdat(pn,:) = interped_data_gd(p,:);   % size  1 x Ncoils (complex)
        else
            Tdat(pn,:) = cal_data(p,:);           % size  1 x Ncoils (complex)
        end

        % if the target signal is too small, we won't use this patch
        if norm(Tdat(pn,:)) < 1e-8
            badPatch(pn) = 1;
        end
        

        % The last layer of the network will be multiplied by the source
        % data (Sdat)
        tmp = cal_data(inds,:)';            % size Nnbrs x Ncoils (complex)
        Sdat(pn,:) = tmp(:)';               % size Nnbrs*Ncoils   (complex)
        
    end
end

% Now reshape the data so it fits into the network for training
% The network looks like this:
% -Input is the 3D coordinates of neighbors :  Npatches x 3*Nnbrs
% -Output is the target signals (complex) : Npatches x Ncoils*2
% -Loss function: compare the output to the target signals.  
% Note that we multiply the last layer (grappa weights) by the Source data 
% this is a complex mulitplication.

Slocs = reshape(nbr_locs, Npatches, 3*Nnbrs);   % neighbor locations 
Sdat = [real(Sdat) imag(Sdat)];                % the source data - used for calculating Loss 
%G_weights = ones(size(Sdat));                  % interpolation weights (still unknown)
Tdat = [real(Tdat) imag(Tdat)];                 % Interpolation targets

% also remove the Infs
[m,n] = find(isinf(Tdat));
badPatch(m)=1;
[m,n] = find(isinf(Sdat));
badPatch(m)=1;

% clean out the NaN's from the training data
% look for Nan's in the first column
badSamps = find(badPatch);

Slocs(badSamps,:) = [];
Sdat(badSamps,:) = [];
Tdat(badSamps,:) = [];
%G_weights(badSamps,:) = [];
Npatches = size(Sdat,1)



% reduce the training set to 1e4 instances
if Npatches > 1e4
    Npatches = 1e4;
    Slocs = Slocs(1:Npatches, :);
    Sdat = Sdat(1:Npatches, :);
    Tdat = Tdat(1:Npatches, :);
    % G_weights = G_weights(1:Npatches, :);
end

Npatches

figure
subplot(311)
hist(Slocs(:), 1000)
title('Neighbors \Delta K xyz' )
subplot(312)
hist(Sdat(:),1000)
title("source signals")
subplot(313)
hist(Tdat(:),1000)
title('Target signals')
drawnow

save training.mat

fprintf('\nStart Training ...\n ');
%%
    

gpuDevice(1)
Lsize = 32
% netParms for perceptron
netParms=[];  % start with a clean variable

netParms.mult1.Weights = (1e-2*dlarray(randn(3*Nnbrs,Lsize)));
%netParms.mult1.Bias = (1e-2*dlarray(randn(1,Lsize)));

netParms.mult2.Weights = (1e-2*dlarray(randn(Lsize,Lsize )));
% netParms.mult2.Bias = (1e-2*dlarray(randn(1,Lsize)));
 
netParms.mult3.Weights = (1e-2*dlarray(randn(Lsize,Lsize)));
% netParms.mult3.Bias = (1e-2*dlarray(randn(1,Lsize)));
 
netParms.mult4.Weights = (1e-2*dlarray(randn(Lsize,Lsize)));
% netParms.mult4.Bias = (1e-2*dlarray(randn(1,Lsize)));
 
netParms.mult5.Weights = (1e-2*dlarray(randn(Lsize,Lsize)));
% netParms.mult5.Bias = (1e-2*dlarray(randn(1,Lsize)));
 
netParms.mult6.Weights = (1e-2*dlarray(randn(Lsize,Lsize)));
% netParms.mult6.Bias = (1e-2*dlarray(randn(1,Lsize)));
% 
netParms.mult7.Weights = (1e-2*dlarray(randn(Lsize,Lsize)));
% netParms.mult7.Bias = (1e-2*dlarray(randn(1,Lsize)));

netParms.mult8.Weights = (1e-2*dlarray(randn(Lsize,Lsize*2)));
%netParms.mult8.Bias = (1e-2*dlarray(randn(1,Lsize*2)));

netParms.mult9.Weights = (1e-2*dlarray(randn(Lsize*2,Lsize*2)));
%netParms.mult9.Bias = (1e-2*dlarray(randn(1,Lsize*2)));

netParms.mult10.Weights = (1e-2*dlarray(randn(Lsize*2,Ncoils*Ncoils*Nnbrs*2)));

%load training.mat

% *****************************************************
% continue training with previous parms ??
% *****************************************************

%load netParms.mat

% *****************************************************

%G_weights = dlarray(G_weights);
Sdat  = gpuArray(Sdat);
Tdat  = gpuArray(Tdat);
Slocs = gpuArray(Slocs);

%
% test the perceptron model function
figure
T_est=model(netParms,Slocs, Sdat);
[loss, gradients]= dlfeval(@modelLoss,netParms, Slocs, Sdat, Tdat);

% figure
% imagesc(squeeze(extractdata(gradients.mult10.Weights)))
% colorbar
% drawnow
% title('Testing: first Gradient of last layer')
%

% define training parameters
iter = 0;
numIter=10;
numIter = 1500;

learnRate = 1e-3;
validationFrequency = 300;


% set up a training monitor object
monitor = trainingProgressMonitor( ...
    Metrics=["TrainingLoss","ValidationLoss"], ...
    Info="Iteration", ...
    XLabel="Iteration");
groupSubPlot(monitor,"Loss",["TrainingLoss","ValidationLoss"]);


% ADAM needs these variables to be present (later)
trailingAvg = [];
trailingAvgSq = [];
minLoss=1e20;

% fractions for training and testing
ftrain = 1:round(Npatches*0.8);
ftest =  round(Npatches*0.8)+1: Npatches;
Ntrain = length(ftrain);

f1= ftrain;
while iter < numIter && ~monitor.Stop
    iter = iter + 1;

    % shuffle the data every few iters.
    % Use only some of the data at each of the shuffles
    %{
    if mod(iter-1,100)==0
        tmp  = randperm(Ntrain);
        f1 = ftrain(tmp(1:500));
    end
    %}

    % Evaluate the model loss and gradients.
    [loss,gradients] = dlfeval(@modelLoss,netParms,Slocs(f1,:), Sdat(f1,:), Tdat(f1,:));
    axis( [-1 1 -1 1]*5 )

    % Update the network netParms using the Adam optimizer.
    [netParms,trailingAvg,trailingAvgSq] = adamupdate(netParms,gradients, ...
        trailingAvg,trailingAvgSq,iter,learnRate);

    if mod(iter-1,20)==0
        % Record the training loss and iter.
        recordMetrics(monitor,iter,TrainingLoss=loss);
        updateInfo(monitor,Iteration=(iter+" of "+numIter));
        monitor.Progress = 100*(iter/numIter);
    end
    % make sure that we hold on to the best network weights
    if (loss < minLoss)
        minLoss=loss;
        savNetParms = netParms;
    end
    
    fprintf('\riter: %d , Loss : %f', iter, loss);
end



netParms=savNetParms;
save netParms.mat netParms
save monitor monitor

Y = model(netParms,Slocs(ftest,:), Sdat(ftest,:));

close

figure


a = gather(Tdat(ftest,:));
b = (extractdata(Y));
rho = corrcoef(a(:), b(:));

plot(a(:),b(:), '.');
xlabel('Truth')
ylabel('Estimate')
axis([-1 1 -1 1 -1 1] *max(abs(a(:))))

title(sprintf('Correlation %f', rho(2,1)))
drawnow
pause(1)

figure

end

