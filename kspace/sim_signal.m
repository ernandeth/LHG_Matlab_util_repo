function kdata = sim_signal( loc , k , object, fov)
% function kdata = sim_signal( xyz_locations , k-trajectory , object_density, fov)
%
% (integrate the signals from all the voxels)
% xyx_locations  - Nx3 matrix
% k-trajectory   - Mx3 matrix
% object         - 3D matrix with densities matching the locations
% fov            - [fovxy fovz]
% returns k-space data gridded into a 3D matrix.  Replaces NaN's with
% zeros.
%
signal = zeros(size(k,1),1);
NPIX = size(loc,1);
fovxy = fov(1);
fovz = fov(2);


for r = 1:NPIX
    %compute the contribution to the signal:
    rvec = loc(r,:);
    fprintf('\rlocation:  %2.2f %2.2f %2.2f -> %2.2f percent   ',rvec(1), rvec(2), rvec(3), r/NPIX*100);
    signal = signal +  object(r)*exp(-i*2*pi*k*rvec' ); 
    %plot(abs(signal));drawnow
    
end
fprintf('\nDone\n');

kmax = max(k(:,2));
kzmax = max(k(:,3));
kzstep = 1/fovz;
kstep = 1/fovxy;

% make a cartesian grid for interpolation of the spiral data
[kxc kyc kzc] = meshgrid(-kmax:kstep: kmax-kstep , ...
     -kmax:kstep: kmax-kstep ,...
     -kzmax:kzstep: kzmax-kzstep);

kdata = griddata3(  k(:,1), k(:,2), k(:,3), signal,kxc, kyc, kzc);
kdata(find(isnan(kdata))) = 0; 

return