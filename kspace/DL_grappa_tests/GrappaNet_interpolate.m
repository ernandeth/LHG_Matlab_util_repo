function [outdata dens] = GrappaNet_interpolate(allCoilNets, signal, ks, dim)
%
% function outdata = grappa_net_interpolate( allCoilNets, signal, ks, dim)
%
% (c) Luis Hernandez-Garcia @University of Michigan 2022
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
%   INPUTS:
%   allCoilNets :  a cell array with the networks.  One network per coil.
%   ks          : [kx, ky, kz] is the trajectory in k-space (cm^-1)
%   signal      : is the signal at each of those ks locations :
%                               (Npts*Nechoes) x Ncoils
%   dim         : is num. pixels in each dimension of the cartesian grid
%                that we are interpolating
%
%   OUTPUT:
%   outdata     : The resulting interpolated cartesian data.
%               ( prod(dim) x Ncoils)
%


%STEP 0:  remove redundant data at center of k-space - (distance ==0)
dist = sqrt(ks(:,1).^2 + ks(:,2).^2 + ks(:,3).^2);
z = find(dist == 0 );
z = z(2:end);
ks(z,:) = [];
signal(z,:) = [];

% STEP 1: set up the cartesian grid and figure out the neighbors
Nnbrs = allCoilNets{1}.Layers(1).InputSize(1);
Ncoils = size(signal,2);
Nsamples = size(signal,1);

fprintf('\nInterpolating FULL cartesian data from non cartesian data using %d neighbors...\n', Nnbrs)

Kmax = max(abs(ks(:)));
dkx = 2*Kmax/dim(1);
Rmax = abs(Kmax/5);
Rmax = dkx*4;   % Radius for interpolation is set to equivalent of 4 kvoxels
% algorithm was trained with 8 kvoxel neighborhood
% but they were half the size.

% The output grid coordinates in k space - units of cm^-1
[xc , yc, zc ]= meshgrid( ...
    linspace(-Kmax,Kmax,dim(1)), ...
    linspace(-Kmax,Kmax,dim(2)), ...
    linspace(-Kmax,Kmax,dim(3)));

Npix = dim(1)*dim(2)*dim(3);  % dimensions of ouput data

outdata = zeros(Npix, Ncoils);  % Cartesian kspace data for all coils
dens = zeros(Npix,1);           % sampling density

% Cartesian kspace data for all coils (for comparison purposes, this is interpolated with griddata)
outdata_grid = zeros(Npix, Ncoils);
%{
for coilnum=1:Ncoils
    fprintf('\n\tInterpolating Cartesian Grid from coil %d (using griddata)',  coilnum)
    % testing for comparison:  use griddata as benchmark
    tmp = griddata(ks(:,1), ks(:,2), ks(:,3), signal(:,coilnum), xc, yc, zc,'linear');
    tmp(isnan(tmp))=eps;
    tmp(isinf(tmp)) = eps;
    outdata_grid(:,coilnum) = tmp(:);
end

fprintf('\nDone with griddata  ... now using Neural Net\n',  coilnum)
%}
%{
% ------------------------------------------------------
% Test:  will this work if I  pass already cartesian data ?
% see if it works for cartesian data as the input
fprintf('\nTest: Replacing spiral data with re-gridded cartesian',  coilnum)
ks = [xc(:) yc(:) zc(:)];
Nsamples = size(ks,1);
signal = outdata_grid;
% --------------------------------------------------
%}
outdata = zeros(Npix, Ncoils);

dist = zeros(Nsamples,1);
nbrlocs = zeros(Nnbrs, 3, 1);
        
for coilnum=1:Ncoils
    
    fprintf('\n\tclearing GPU...')
    reset(gpuDevice(1));
    wait(gpuDevice(1))

    fprintf('\n\tInterpolating Cartesian Grid from coil %d (using Neural Net)',  coilnum)
    mynet = allCoilNets{coilnum};
    
    % for each of the targets on cartesian grid (all locations)
    for p =  1:Npix
        
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
        % first find distances from target voxels (cartesian gris) to all sample locs
        % and then consider only those that are close enough (<Rmax)

        
        for n=1:Nsamples
            dist(n) = sqrt(...
                (xc(p) - ks(n,1))^2 + ...
                (yc(p) - ks(n,2))^2 + ...
                (zc(p) - ks(n,3))^2);
        end
        
        inds = find(dist < Rmax);
        % save the sampling density for later use
        dens(p) = length(inds);
        
        if length(inds) < Nnbrs
            % if no neighbors are close enough, we put in a zero
            W =0;
            outdata(p,coilnum) = 0;
            %fprintf('\rNo neighbors at %f %f %f    :( ', xc(p), yc(p), zc(p));
        else
            % Choose a random set of Nnbrs to make an interpolation 'patch'
            %tmp = randi(length(inds), Nnbrs,1);
            %inds = inds(tmp);
            
            % Choose the closest points to make an interpolation patch
            [dist inds] = sort(dist);
            inds = inds(1:Nnbrs);
            
            patch_data = signal(inds,:); % size Nnbrs x Ncoils (complex)
            %patch_data = patch_data';  %   <---  Test ????
            patch = patch_data(:);
            
            % Now calculate locations of neighbors RELATIVE to the target
            % for this patch
            nbrlocs(:,:,1) = [ ks(inds,1)-xc(p) ,  ks(inds,2)-yc(p) ,   ks(inds,3)-zc(p)];
            
            % STEP 3: compute the interpolation weights using the network
            % Use the neural net to produce the corresponding COMPLEX weights:
            % (the output of the network has all the reals followed by  the
            % imaginaries - must merge into complex numbers)
            
            W = predict(mynet , nbrlocs);
            W = complex(W(1:end/2), W(end/2+1:end));
            
            % W = W/norm(W);  % should I normalize the weights?
            
            % STEP 4: Interpolate the signals at the desired Cartesian location
            % as a weighted sum of the neighbors at each coil
            outdata(p,coilnum) = W*patch;
            
            
            %
            if mod(p,1e3)==0
                fprintf('\r\t\tCoil %d Kvoxel %d of %d', coilnum, p, Npix);
                %{
                subplot(221)
                hold off
                plot3(ks(inds,1), ks(inds,2), ks(inds,3), 'ro' );
                hold on
                plot3(xc(p), yc(p), zc(p), 'kx' );
                axis ([-1 1 -1 1 -1 1]*Kmax)
                axis square
                title('Neighbors')
                
                subplot(222)
                hold off; plot(abs(W),'r')
                hold on; plot(abs(patch/max(patch)),'b')
                legend('weights', 'signal');
                
                subplot(223)
                lightbox(reshape(dens, dim));
                title('Density')
               
                whos W patch
                toc
            
                drawnow
                %}
             end
            %}
        end
    end  % kvoxel loop
                
    clear mynet

    %}
    %subplot(224)
    %
    % compare the result of my interpolation to the griddata interpolation
    tmp = outdata(:,coilnum);
    tmp = reshape(tmp, dim);
    
    tmp2 = outdata_grid(:,coilnum);
    tmp2 = reshape(tmp2, dim);
   
    rho = corrcoef(abs(tmp(:)), abs(tmp2(:)));
    %{
    figure
    subplot(221)
    plot(abs(tmp(:)), abs(tmp2(:)),'.');
     fprintf('correlation beteen linear and DL interp, rho = %f' rho(1,2)));
    
    subplot(223)
    lightbox(log(abs(tmp)));
    title('NN interpolation ')
    
    subplot(224)
    lightbox(log(abs(tmp2)));
    title('Gridding interpolation')
    
    drawnow
    %}
    clear tmp tmp2
    title(sprintf('correlation between griddata and DLinterp for coil %d, rho = %f', coilnum, rho(1,2)));
end  % coil loop

return