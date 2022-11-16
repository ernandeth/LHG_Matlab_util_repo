function [outdata nb_save wt_save] = myinterp3d(indata, inx, iny, inz, outx, outy, outz, N_nbors)
% function [outdata nb_weights] = myinterp3d(indata, inx, iny, inz, outx, outy, outz, N_nbors)

outloc = [ outx(:) outy(:) outz(:)];
Noutlocs  = size(outloc,1);
outdata = zeros(Noutlocs,1);

% separation in output k-space
deltak = abs(outx(1,1,1) -outx(1, 2,1)); 

inloc = [inx(:) iny(:) inz(:)];

% % adjusting for acquisition delays
% inloc = inloc(1 : end-26-10,:);
% indata = indata(27 : end-10);
% 
% % we can probably do it with just half the data
% inloc = inloc(1:2:end, :);
% indata = indata(1:2:end);
% %%%%%%%%%%%
Ninlocs = size(inloc,1);

Nweights = (1+N_nbors)^3 ;
weights = zeros(1, Nweights);
nbors = zeros(1, Nweights);
dists = zeros(1, Nweights);

% we will store the position of the neighbors and the weights for 
% interpolating subsequent images
nb_save = zeros(Noutlocs, Nweights);
wt_save = zeros(Noutlocs, Nweights);

alldists = zeros(1, Ninlocs);
here = zeros(Ninlocs,3);

% the normalization constant for the gaussian weights =  area of the gaussian
Gn = (1/N_nbors*deltak*sqrt(2*pi));

for n=1:Noutlocs
    
    if ~mod(n,10000)
        fprintf('\rworking on kloc %d of %d ...', n, Noutlocs);
    end
    
    nbors = zeros(1, Nweights);
    dists = zeros(1, Nweights);

    % first calculate the distance to all the existing data points
    % and choose the closest ones
    tmp = outloc(n,:);
    here = repmat(tmp, Ninlocs,1);
    
    alldists = sqrt(sum((inloc - here).^2 , 2) ) ;
    
   
    % find the neighbors (anything closer than N_nbors)
    inds = find(alldists <= N_nbors*deltak);    
    tmp_dists = alldists(inds);
    Nfound = length(inds);
    
    % if I have too many, just keep the colsest ones for the interpolation
    % (sort that neighborhood and keep the closest Nweights)
    if Nfound > Nweights
            [tmp_dists p] = sort(tmp_dists);
            
            inds = inds(p(1:Nweights));
            tmp_dists = tmp_dists(1:Nweights);
            Nfound = Nweights;
    end
        
    
    % do the interpolation
    if length(inds) < N_nbors/2
        % if I don't have enough neighbors for a decent average
        % just put in a zero and move on
        weights = zeros(1,Nweights);
        nbors = zeros(1,Nweights);
        outdata(n) = 0;
       
    else
        dists(1:Nfound) = tmp_dists';
        nbors(1:Nfound) = inds' ;
        
        % now use them for
        % the weighting kernel calculation
        weights = sinc(dists/deltak);
        %weights = Gn *exp( - (dists / (2*N_nbors*deltak)).^2 );
        
        % the voxel's signal is the weighted average of the neighbors
        outdata(n) = sum(indata(nbors(1:Nfound) ) .* weights(1:Nfound)) / Nfound ;

        % store that info for the next frame
        nb_save(n, : )= nbors;
        wt_save(n, : ) = weights;
%         

    end
end

save nbrs_weights.mat  nb_save wt_save

% figure
% outdata = reshape(outdata, length(outx), length(outy), length(outz));
% lightbox(abs(outdata));

return


% example code:
ndata = 200;
lin = linspace(-10,10,20);
[outx, outy, outz] = meshgrid(lin, lin, lin);

r = linspace(0,10,ndata);
theta = linspace(0,10*pi, ndata);
inx = real(r.*exp(i*theta));
iny = imag(r.*exp(i*theta));
inz = 10*sin(linspace(0,pi/2, ndata));
plot3(inx, iny, inz);

indata = exp(-linspace(0,1,ndata));
%indata = ones(1,ndata);
N_nbors = 2;
[outdata nb_weights] = myinterp3d(indata, inx, iny, inz, outx, outy, outz, N_nbors);
