function [im,fov,tr] = nii2matrix(name)
% This function is a wrapper of read_nii_img, that can read an image from a
% file and automatically reshape/descale the data based on info from the
% header

    [raw,h] = read_nii_img(name);
    
    % Get dimensions:
    tr = h.pixdim(5);
    dim = h.dim(2:4);
    fov = h.dim(2:4).*h.pixdim(2:4);
    nframes = h.dim(5);
    
    % Reshape data into matrix with time on 4th dimension:
    im = permute(reshape(raw,[nframes,dim]),[2:4,1]);
    
    % Descale data:
    if ~isempty(h.scl_slope), im = im*h.scl_slope; end
    if ~isempty(h.scl_inter), im = im + h.scl_inter; end
    
end

