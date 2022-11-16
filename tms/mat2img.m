Eabs_nometa = sqrt(real(Ebrain(:,:,:,1)).^2 + ...
    real(Ebrain(:,:,:,2)).^2 + real(Ebrain(:,:,:,3)).^2); 

E = Eabs_nometa;

hdr = define_avw_hdr;
hdr.sizeof_hdr = 348;
hdr.dims = 3;
hdr.tdim=1;
hdr.xdim = 80;
hdr.ydim = 87;
hdr.zdim = 44;
hdr.xsize = 0.018;
hdr.ysize = 0.018;
hdr.zsize = 0.03;
hdr.datatype = 16;
hdr.bits = 32;
 write_hdr('mysigma.hdr', hdr);
 write_img('mysigma.img', sigma, hdr);
% 
% 
% 
% Efield=sqrt(E(:,:,:,1).^2+E(:,:,:,2).^2+E(:,:,:,3).^2);
write_hdr('myEfield_nometa.hdr', hdr);
write_img('myEfield_nometa.img', E, hdr);