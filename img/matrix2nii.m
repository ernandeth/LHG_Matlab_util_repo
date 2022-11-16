function matrix2nii(name,data,fov,tr,doScl)
% This function is a wrapper of write_nii, that takes in a 4D matrix and
% writes it to file in signed 16-bit int format based on necessary
% parameters. It can utilize the full dynamic range of the 16-bit by
% passing in doScl = 1

    if nargin<6 || isempty(doScl)
        doScl = 1;
    end

    if nargin<5 || isempty(tr)
        tr = 1;
    end

    % Signed 16-bit integer format scaling factor:
    DACFactor = 2^15-1;
    
    % get dimensions
    dim = [size(data,1),size(data,2),size(data,3)];
    nframes = size(data,4);
    
    % define header
    h = define_nii_hdr();
    h.dim = [4 dim nframes 0 0 0];
    h.pixdim = [4 fov./dim tr 0 0 0];
    h.datatype = 4;
    h.bitsperpixel = 16;
    
    % scale to full dynamic range and write scaling value to file
    if doScl
        y = data;
        y_min = min(y,[],'all'); y_max = max(y,[],'all');
        m = 2*DACFactor/(y_max - y_min);
        x = m*y - (m*y_min + DACFactor);
        h.scl_inter = (m*y_min + DACFactor)/m;
        h.scl_slope = 1/m;
        data = x;
    end
    
    % write out to file
    write_nii(name,data(:),h,0);
    
end

