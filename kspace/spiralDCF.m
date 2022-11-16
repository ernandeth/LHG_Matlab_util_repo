function kdens = spiralDCF(ks,N,R)
%function SPIRALDCF(ks,N,R)
%|
%|  Function that calculates the sampling density compensation function
%|      based on k-space trajectory by averaging distance to neighboring
%|      samples
%|
%|  Inputs:
%|      ks: kspace trajectory [ks_x ks_y ks_z]
%|      N: number of neighbors in density calcualtion neighborhood
%|      R: exponent on distance function
%|

    f = @(x) x.^R;
    fprintf('Creating density compensation...\n')
    kdens = zeros(size(ks,1),1);
    parfor isample = 1:size(ks,1)
        ksub = ks - ks(isample,:);
        dists = vecnorm(ksub,2,2);
        kdens(isample) = mean(mink(dists,N));
    end
    kdens = (kdens - min(kdens,[],'all'))./(max(kdens,[],'all') - min(kdens,[],'all'));
    kdens = f(kdens);
    kdens(isinf(kdens)) = 1; kdens(isnan(kdens)) = 0;

end

