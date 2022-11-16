dim = ones(3,1)*65;
Kmax = 65/24/2;
Ncoils = 8;
Nnbrs = 21;

% The output grid coordinates in k space
[xc , yc, zc ]= meshgrid( ...
    linspace(-Kmax,Kmax,dim(1)), ...
    linspace(-Kmax,Kmax,dim(2)), ...
    linspace(-Kmax,Kmax,dim(3)));

[xi, yi, zi] = meshgrid([1:dim(1)], [1:dim(2)], [1:dim(3)]);

weights_grid = zeros(dim(1), dim(2), dim(3), Ncoils);

Npix = dim(1)*dim(2)*dim(3);
    
load GRAPPAnet.mat

for p = 1:10000
    inds = randi(65, Nnbrs, 3);
    nbrlocs = Kmax*(inds-33)/32;
    inbrlocs =inds;
    
    mynet = allCoilNets{1};
    W = predict(mynet , nbrlocs(:,:,:));
    W = complex(W(1:end/2), W(end/2+1:end));
    W = reshape(W, Nnbrs, Ncoils);
    
    for n=1:Nnbrs
        weights_grid(inbrlocs(n,1), inbrlocs(n,2), inbrlocs(n,3)) = W(n,1);
    end
end

lightbox(log(abs(weights_grid(:,:,:,1))));