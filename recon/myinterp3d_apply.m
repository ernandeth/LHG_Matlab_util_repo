function [outdata ] = myinterp3d_apply(indata,  outx, outy, outz, nbrs, wts)
% function [outdata nb_weights] = myinterp3d(indata,  outx, outy, outz, N_nbors)

outloc = [ outx(:) outy(:) outz(:)];
Noutlocs  = size(outloc,1);
outdata = zeros(Noutlocs,1);



for n=1:Noutlocs;
    
    if ~mod(n,10000)
        fprintf('\rapplying interpolation on kloc %d of %d ...', n, Noutlocs);
    end
    
    weights = wts(n,:);
    neighbors = nbrs(n,:);
    % remove the zero neighbors (ie - no neighbors)
    inds = find(neighbors);
    Nweights = length(inds);
    
    % the voxel's signal is the weighted average of the neighbors
    if Nweights==0
        outdata(n) = 0;
    else
        
        outdata(n) = sum(indata(neighbors(inds)) .* weights(inds)) / Nweights ;
    end
end

outdata = reshape(outdata, size(outx,1), size(outy,2), size(outz,3));

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
indata = ones(1,ndata);
[outdata nb wt] = myinterp3d(indata, inx, iny, inz, outx, outy, outz, 2);
figure
outdata = myinterp3d_apply(indata,  outx, outy, outz, nb, wt)

