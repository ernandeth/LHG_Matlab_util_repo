function [interped_data interped_data_gd dens] = GrappaNet_interpolate(modelParms, signal, ks, dim)
%
% function [interped_data interped_data_gd dens ]= grappa_net_interpolate( allCoilNets, signal, ks, dim)
%
% (c) Luis Hernandez-Garcia @University of Michigan 2022
%
%
% Interpolate cartesian data (targets) from
% non-cartesian data (neighbors) using a neural net.
%
% Target(r, c) = sum{ neighbors(delta_r, c) * weight(delta_r, c)}
%                   over all delta_r and coils
%
% The Interpolation Weights are calculated by a neural net
% based on the positions of the neighbors relative to the targets.
%
% We assume that there is a function that determines the value of interpolation
% kernel (the weights) at the locations of the neighbors in 3D.
%
% This function was learned by a neural network (different code).
% It was trained from fully sampled data in a calibration region.
% The network's inputs are the relative locations
% of the interpolation neighbors, and their signal values
% the output is the signal at the target location in kspace for each coil.
%
%   INPUTS:
%   allCoilNets :  a cell array with the networks.  One network FOR ALL coils.
%   ks          : [kx, ky, kz] is the trajectory in k-space (cm^-1)
%   signal      : is the signal at each of those ks locations :
%                               (Npts*Nechoes) x Ncoils
%   dim         : is num. pixels in each dimension of the cartesian grid
%                that we are interpolating
%
%   OUTPUT:
%   interped_data     : The resulting interpolated cartesian data.
%               ( prod(dim) x Ncoils)
%


%STEP 0:  remove redundant data at center of k-space - (distance ==0)
% keep only one of those points
R = sqrt(ks(:,1).^2 + ks(:,2).^2 + ks(:,3).^2);
z = find(R == 0 );
z = z(2:end);
ks(z,:) = [];
signal(z,:) = [];

% STEP 1: set up the cartesian grid and figure out the num. of neighbors
Nnbrs = size(modelParms.mult1.Weights,1)/3;
Ncoils = size(signal,2);
Nsamples = size(signal,1);

fprintf('\nInterpolate Cartesian data from non-cartesian data using %d neighbors...\n', Nnbrs)

Kmax = max(R(:));
dkx = 2*Kmax/dim(1);
Rinterp = 3*dkx;   % Radius for interpolation is set to equivalent of 3 kvoxels


%%%%% reduce the size of the output k-space data set to speed up testing %%%%
% Kmax = Kmax/2;
% dim = dim/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The output grid coordinates in k space - units of cm^-1
[xc , yc, zc ]= meshgrid( ...
    linspace(-Kmax,Kmax,dim(1)), ...
    linspace(-Kmax,Kmax,dim(2)), ...
    linspace(-Kmax,Kmax,dim(3)));

Npix = dim(1)*dim(2)*dim(3);  % dimensions of output data
kc = ([xc(:) yc(:) zc(:)]);   % the coordinates of the cartesian grid as an Nx3

% Cartesian kspace data for all coils - output of the neural net estimation
interped_data = zeros(Npix, Ncoils);  
dens = zeros(Npix,1);                 % sampling density
% Cartesian kspace data for all coils (for comparison purposes, this is interpolated with griddata)
interped_data_gd = zeros(Npix, Ncoils);

fprintf('\nTry griddata for reference ...\n')

parfor coilnum=1:Ncoils
    fprintf('\nInterpolating Cartesian Grid from coil %d using griddata',  coilnum)
    % testing for comparison:  use griddata as benchmark
    tmp = griddata(ks(:,1), ks(:,2), ks(:,3), signal(:,coilnum), xc, yc, zc,'natural');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    interped_data_gd(:,coilnum) = tmp(:);
end

fprintf('\n...Done with griddata. Now using Neural Net ...\n')


%}
%{
% ------------------------------------------------------
% Test:  will this work if I  pass already cartesian data ?
% see if it works for cartesian data as the input
fprintf('\nTest: Replacing spiral data with re-gridded cartesian',  coilnum)
ks = [xc(:) yc(:) zc(:)];
Nsamples = size(ks,1);
signal = interped_data_gd;
% --------------------------------------------------
%}

%{
fprintf('\n\tclearing GPU...')
reset(gpuDevice(1));
wait(gpuDevice(1))
fprintf('GPU cleared')
%}
randomNbrs = 1;

% Allocate space for the network inputs and outputs:

% source signals
Sdat = (zeros(Npix, Nnbrs*Ncoils));   
% source positions relative to the target
Slocs = (zeros(Npix, Nnbrs,3));
% OUTPUT of the network: target signals (interpolated)  [real imag] 
Tdat = (zeros(Npix, Ncoils*2));

% buffer for the distances from target to neightbors
dist = (zeros(Nsamples,1));
% buffer for the neighbor locations  ... same as Slocs
nbr_locs = (zeros(Nnbrs, 3, 1));

fprintf('\nBegin collecting Neighbor Signals and locations for interpolation...\n')


% for each of the targets on cartesian grid (all locations)
for p = 1:Npix

    % identify neighbors for interpolation :
    % first find distances from target voxels (cartesian grid) to all
    % (noncartesian) sample locs
    dist = sqrt(sum( (kc(p,:)-ks).^2 ,2));
    % and then consider only those that are close enough (<Rinterp)
    nbr_inds = find(dist < Rinterp  & dist > 0);

    if length(nbr_inds) >= Nnbrs

        if randomNbrs==1
            % Choose a random set of Nnbrs to make an interpolation 'nbr_signals'
            tmp = randperm(length(nbr_inds), Nnbrs);
            inds = nbr_inds(tmp);

        else
            % ... OR ... Choose the closest points to make an interpolation nbr_signals
            [dist inds] = sort(dist);
            inds = inds(1:Nnbrs);
        end

        % Now calculate location vectors of neighbors RELATIVE to the target
        % for this patch of data
        nbr_locs = ks(inds,:) - kc(p,:);
                
        % non-cartesian data corresponding to neighbors
        nbr_signals = signal(inds,:)'; % size Nnbrs x Ncoils (complex)
                
        % sampling density: how close are these neighbors?
        dens(p) = 1/mean(vecnorm(nbr_locs));

        % We now have a patch for interpolation: 
        % neighbor signals and their locations 
        % Store them into the input matrices for the Neural net
        Sdat(p,:) = nbr_signals(:)';
        Slocs(p,:,:) = nbr_locs;

        %{
        % Show locations of neightbors and targets for debugging
        if mod(p,200)==0
            fprintf('\rGetting Neighbors for Kvoxel %d of %d',  p, Npix);
            
            subplot(221)
            hold off
            plot3(ks(inds,1), ks(inds,2), ks(inds,3), 'r.' );
            hold on
            plot3(xc(p), yc(p), zc(p), 'k.' );
            axis ([-1 1 -1 1 -1 1]*Kmax)
            axis square
            title('Neighbors')

            subplot(222)
            %hold off; plot(abs(Sdat),'r')
            imagesc(abs(nbr_signals./max(nbr_signals)))
            ylabel('coils')
            xlabel('neighbors')
            colormap gray; colorbar
            title('neighbor signals')

            drawnow
            
        end
        %}
        
    end % do calculations if there are enough neighbors
end  % kvoxel loop

fprintf('\n...Finished collecting Neighbors \n')

Sdat = [real(Sdat) imag(Sdat)];
Slocs = reshape(Slocs, Npix, Nnbrs*3);

% convert to dlarray in order to use DL toolbox
% maybe gpuArray?
Sdat = dlarray(Sdat);
Slocs = dlarray(Slocs);
Tdat = dlarray(Tdat);

fprintf('\nBegin interpolation using Network ...')

% Compute the target signals using the network
% The input is the nbr_locs and the corresponding nbr_signals as
% the additional parameter
% The network updates internal weghts, such that the last layer
% contains the apporpriate interpolation weights.  These are
% multiplied by the corresponsing signals within the network.
% The output of the network is the target location's signals
% in all coils.
% (the output of the network has all the reals followed by  the
% imaginaries - must merge into complex numbers)
bufsize=10000;
for p=1:bufsize:Npix
    inds = p:min(Npix, p+bufsize-1);
    fprintf('\nInterpolating pixel %d of %d', p, Npix);
    Tdat(inds,:) = model(modelParms, Slocs(inds,:), Sdat(inds,:));
end
fprintf('\n...Finished interpolation\n')

%{
Rmax = sqrt(sum(kc.^2,2 ));
p=find(Rmax>Kmax);
Tdat(p,: )=0;
%}

%
fprintf('\nSecond Pass to fill in holes ...')
% Second Pass:  fill the gaps
for p = 1:Npix

    dist = sqrt(sum( (kc(p,:)-ks).^2 ,2));
%
    % See if any of the points coincided with the cartesian grid
    % ie - no need to interpolate:
    nbr_inds = find(dist < 0.01*dkx);
    if ~isempty(nbr_inds)
        tmp = mean(signal(nbr_inds, :), 1);
        Tdat(p,:) = [real(tmp) imag(tmp)]; 
        if ~isempty(find(isnan(tmp)))
            fprintf('\n Overlap @ pixel %d - average neighbors', p)
        end
    end
%
    % See which points didn't have enough neighbors.  We will use a simple
    % average of the neighbors for interpolation in these cases
    nbr_inds = find(dist < Rinterp  & dist>0);
    if isempty(nbr_inds)
            Tdat(p,:) = 0;

    else if length(nbr_inds) < Nnbrs
            tmp = mean(signal(nbr_inds, :), 1);
            Tdat(p,:) = [real(tmp) imag(tmp)];

            if ~isempty(find(isnan(tmp)))
                fprintf('\n Too few nbrs at pixel %d - average the neighbors we have', p)
            end
    end
    end

end
fprintf('\n ... Second Pass complete.')
%}

% interped_data is a dlarray - must extract the data from it first
interped_data = extractdata(complex(Tdat(1:end/2), Tdat(end/2+1:end)));
interped_data = (reshape(interped_data, Npix, Ncoils));

% show k-space data of one of the coils to see if it works 
figure
subplot(211)
lightbox(reshape(log(dens), dim));
title('Density ()')

subplot(212)
lightbox(reshape(log(abs(interped_data(:,5))), dim));
title('coil 5 Kspace data (log)')



% compare the result of my interpolation to the griddata interpolation
for coilnum=1:Ncoils

    tmp = interped_data(:,coilnum);
    tmp = reshape(tmp, dim);

    tmp2 = interped_data_gd(:,coilnum);
    tmp2 = reshape(tmp2, dim);

    rho = corrcoef(abs(tmp(:)), abs(tmp2(:)))
    %{
    figure
    subplot(221)
    plot(abs(tmp(:)), abs(tmp2(:)),'.');
    title(sprintf('\nch %d Corr. grid v. DL = %f' , coilnum, rho(1,2)))

    subplot(223)
    lightbox(log(abs(tmp)));
    title('NN interpolation ')

    subplot(224)
    lightbox(log(abs(tmp2)));
    title('Gridding interpolation')

    drawnow
    %}
    clear tmp tmp2
end  % coil loop

return