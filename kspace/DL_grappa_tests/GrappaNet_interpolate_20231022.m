function [interped_data interped_data_gd dens] = GrappaNet_interpolate(allCoilNets, signal, ks, dim)
%
% function [interped_data interped_data_gd dens ]= grappa_net_interpolate( allCoilNets, signal, ks, dim)
%
% (c) Luis Hernandez-Garcia @University of Michigan 2022
%
% This version tries to compute all the coils together using a single
% network (2023.10.22)
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
% of the interpolation neighbors, the output is the weight for each
% neighbor and each coil. 
% 
% 2023.01.02 The training is done in non-cartesian data
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
dist = sqrt(ks(:,1).^2 + ks(:,2).^2 + ks(:,3).^2);
z = find(dist == 0 );
z = z(2:end);
ks(z,:) = [];
signal(z,:) = [];

% STEP 1: set up the cartesian grid and figure out the neighbors
Nnbrs = allCoilNets.Layers(1).InputSize(1);
Ncoils = size(signal,2);
Nsamples = size(signal,1);

fprintf('\nInterpolate Cartesian data from non-cartesian data using %d neighbors...\n', Nnbrs)

Kmax = max(abs(ks(:)));
dkx = 2*Kmax/dim(1);
%Rmax = abs(Kmax/5);
Rmax = dkx*4;   % Radius for interpolation is set to equivalent of 4 kvoxels
% algorithm was trained with 8 kvoxel neighborhood
% but they were half the size.

%%%%% reduce the size of the k-space data set for testing %%%%
Kmax = Kmax/4;
dim = dim/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The output grid coordinates in k space - units of cm^-1
[xc , yc, zc ]= meshgrid( ...
    linspace(-Kmax,Kmax,dim(1)), ...
    linspace(-Kmax,Kmax,dim(2)), ...
    linspace(-Kmax,Kmax,dim(3)));

Npix = dim(1)*dim(2)*dim(3);  % dimensions of ouput data

interped_data = zeros(Npix, Ncoils);  % Cartesian kspace data for all coils
dens = zeros(Npix,1);                 % sampling density

% Cartesian kspace data for all coils (for comparison purposes, this is interpolated with griddata)
interped_data_gd = zeros(Npix, Ncoils);
%
fprintf('\nTry griddata first ...\n')

parfor coilnum=1:Ncoils
    fprintf('\n\tFor Reference:  Interpolating Cartesian Grid from coil %d using griddata',  coilnum)
    % testing for comparison:  use griddata as benchmark
    tmp = griddata(ks(:,1), ks(:,2), ks(:,3), signal(:,coilnum), xc, yc, zc,'linear');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    interped_data_gd(:,coilnum) = tmp(:);
end

fprintf('\nDone with griddata  ... now using Neural Net ...\n')
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
interped_data = zeros(Npix, Ncoils);

dist = zeros(Nsamples,1);
nbrlocs = zeros(Nnbrs, 3, 1);



fprintf('\n\tclearing GPU...')
reset(gpuDevice(1));
wait(gpuDevice(1))
fprintf('GPU cleared')



% for each of the targets on cartesian grid (all locations)
for p = 1:Npix % dim(1)*dim(2)*dim(3)/2 + dim(1)*dim(2)/2 + dim(1)/2   % 1:Npix

    %         % ------------------ Test   ---------------------
    %         % interpolate only the calibration region (this
    %         is where the networks were trained
    %         if ((abs(xc(p)) < Kmax/2 ) & ...
    %             (abs(yc(p)) < Kmax/2 ) & ...
    %             (abs(zc(p)) < Kmax/2 ))

    %         if mod(p,1e4)==0
    %             tic
    %         end

    % STEP 2 - identify neighbors for interpolation :
    % first find distances from target voxels (cartesian grid) to all 
    % (noncartesian) sample locs
    % and then consider only those that are close enough (<Rmax)
    parfor n=1:Nsamples
        dist(n) = sqrt(...
            (xc(p) - ks(n,1))^2 + ...
            (yc(p) - ks(n,2))^2 + ...
            (zc(p) - ks(n,3))^2);
    end

    nbr_inds = find(dist < Rmax  & dist>0);
    % save the sampling density for later use
    dens(p) = length(nbr_inds);

    if length(nbr_inds) < Nnbrs
        % if no neighbors are close enough, we put in a zero
        W =0;
        interped_data(p,:) = 0;
        %fprintf('\rNo neighbors at %f %f %f    :( ', xc(p), yc(p), zc(p));
    else
        % Choose a random set of Nnbrs to make an interpolation 'patch'
%         tmp = randi(length(nbr_inds), Nnbrs,1);
%         inds = nbr_inds(tmp);

        % ... OR ... Choose the closest points to make an interpolation patch
        [dist inds] = sort(dist);
        inds = inds(1:Nnbrs);

        patch_data= signal(inds,:); % size Nnbrs x Ncoils (complex)
        patch = patch_data(:)';

        % Now calculate locations of neighbors RELATIVE to the target
        % for this patch
        nbrlocs(:,:,1) = [ ks(inds,1)-xc(p) ,  ks(inds,2)-yc(p) ,   ks(inds,3)-zc(p)];

        %dens(p) = mean(vecnorm(nbrlocs,2,2));

        % STEP 3: compute the interpolation weights using the network
        % Use the neural net to produce the corresponding COMPLEX weights:
        % (the output of the network has all the reals followed by  the
        % imaginaries - must merge into complex numbers)

        % choose the right net for each channel and use it to calculate
        % the interp. weights
        mynet = allCoilNets;

        %W = predict(mynet , nbrlocs);

        gpuNbrlocs = gpuArray(nbrlocs);

        W = predict(mynet , gpuNbrlocs);
        W = complex(W(1:end/2), W(end/2+1:end));

        % reshape the Weights into the original size:
        % [ Nnbrs*Ncoils  x  Ncoils]
        W = reshape(W, [ Nnbrs*Ncoils ,  Ncoils]);

        % debugging - how does it look if they have equal weight?
        %W = ones(1,length(patch));

        % W = W/norm(W);  % should I normalize the weights? - NO .
        % Makes things worse

        % STEP 4: Interpolate the signals at the desired Cartesian location
        % as a weighted sum of the neighbors at each coil
        interped_data(p,:) = patch*W;

        %
        if mod(p,1e3)==0
            fprintf('\r\t\tKvoxel %d of %d',  p, Npix);
            %
            subplot(221)
            hold off
            plot3(ks(inds,1), ks(inds,2), ks(inds,3), 'r.' );
            hold on
            plot3(xc(p), yc(p), zc(p), 'ko' );
            axis ([-1 1 -1 1 -1 1]*Kmax*2)
            axis square
            title('Neighbors')

            subplot(222)
            hold off; plot(abs(W),'r')
            hold on; plot(abs(patch/max(patch)),'b')
            legend('weights', 'signal');

            subplot(223)
            lightbox(reshape((dens), dim));
            title('Density ()')

            subplot(224)
            lightbox(reshape(log(abs(interped_data(:,5))), dim));
            title('coil 5 Kspace data (log)')

            % whos W patch
            nbrlocs

            drawnow
            %}

            %}
        end % Every 1e3 voxels : display results
    end % do calculations if there are enough neighbors
end  % kvoxel loop

clear mynet

%}
%subplot(224)
%

% compare the result of my interpolation to the griddata interpolation
for coilnum=1:Ncoils

    tmp = interped_data(:,coilnum);
    tmp = reshape(tmp, dim);

    tmp2 = interped_data_gd(:,coilnum);
    tmp2 = reshape(tmp2, dim);

    rho = corrcoef(abs(tmp(:)), abs(tmp2(:)));
    %
    figure
    subplot(221)
    plot(abs(tmp(:)), abs(tmp2(:)),'.');
    fprintf('correlation beteen linear and DL interp, rho = %f' , rho(1,2));

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