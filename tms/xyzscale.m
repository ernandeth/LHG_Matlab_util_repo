function result = xyzscale()
% function result = xyzscale()
% converts data from pixels to cm units

[n p]=uigetfile('*.hdr', 'Select image template');
h=read_hdr(strcat(p,n));
scale(1)=h.xsize/10;
scale(2)=h.ysize/10;
scale(3)=h.zsize/10;


[n p]=uigetfile('*','Select data to rescale');
name=strcat(p,n);
data=read_mat(name,4)

result(:,1) = data(:,1) * scale(1);
result(:,2) = data(:,2) * scale(2);
result(:,3) = data(:,3) * scale(3);
result(:,4) = data(:,4);
result
name=strcat('cm-',n)
write_mat(result,strcat(p,name));

return






