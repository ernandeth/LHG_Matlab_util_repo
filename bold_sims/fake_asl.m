
h=read_hdr('perf004.hdr');
im=read_img2(h,'perf004.img');
show(im(:,:,3));


m=zeros(100,1);
%m(5:8)=.8;
%m(15:18) = 0.8;
%m(30:33) = 0.8
m(1:33)= spm_hrf(1);
m(31:63) =spm_hrf(1);
m(61:93) =spm_hrf(1);
%m = m+1;

m(2:2:end)=- 0.5* m(2:2:end)-0.1; 
plot(m)

for i=1:size(m,1)
	
      im2 = im;
      im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, h.zdim/2:h.zdim) = im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75,h.zdim/2:h.zdim ) * (1 + m(i)* 0.1 ) ;
      str=sprintf('asl_%04d.img',i);
      %im2 = im2 *(1 +  mod(i,2)*0.1);
      write_img_data(str,im2,h);
      
      str=sprintf('asl_%04d.hdr',i);
      write_hdr(str,h);
      
      
end

return

datasize=15
x=30; y=30; z=3;

index= h.xdim*h.ydim* (z-1) + h.xdim*(y-1) + x-1 ;
for i=1:datasize
   str=sprintf('asl_%04d.img',i)
   pFile = fopen(str,'rb');
   fseek(pFile,index*2,'bof');
   series(i)=fread(pFile,1,'int16');
   fclose(pFile);
end

