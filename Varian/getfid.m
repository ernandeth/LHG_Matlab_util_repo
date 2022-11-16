function [kspace_data,image_data,image_abs,image_phase] = getfid()

[file_name,file_path]=uigetfile('*.*','Choose Raw FID file');
cd(file_path);
filename=strcat(file_path,file_name);

[kspace,np,nt] = ReadVarian2D(filename);
kspace_data = fliplr(flipud(reshape(kspace,np/2,nt)));


image_data = ifftshift(ifft2(kspace_data));
image_abs=abs(image_data);
image_phase=angle(image_data);


return
