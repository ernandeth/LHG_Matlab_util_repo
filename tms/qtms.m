%The following script is to reconstruct magnitude images and field maps using gsp21a. If any program fails, ask Luis about setting up paths properly. The necessary programs are : Luis' m file collection, SPM2 package, Sangwoo's m file collection


!rm -f *.img
!rm -f *.hdr
!rm -f ref*
!rm -f *.mat
!rm -f vol*
!rm -f rvol*

!~/bin/gsp21a -m -l -fx -fy  Popen1
!~hernan/scripts/ref-stack open1 128 128 60 1 1 1 ref.128.s

!~/bin/gsp21a -m -l -fx -fy  Pon1
!~hernan/scripts/ref-stack on1 128 128 60 1 1 1 ref.128.s
                                                                                
!~/bin/gsp21a -A -l -fx -fy -n 128 Pon1 Popen1
                                                                                

v=dir('vol*.img');
anat = read_img2(v(1).name);
mask=zeros(size(anat));
threshold = 0.1*max(anat(:));
mask(find(abs(anat) > threshold )) = 1;

% this part is to align the field map orientation to the magnitude images.
open1=read_img2('open1.img');
on1=read_img2('on1.img');
hdr1 = read_hdr('on1.hdr');

% calculate the induced Bz field from the phase difference map. 
mapdel=2.5e-3; % mapdelay in sec
gamma_bar=4258; % Hz / Gauss
Bz_1= mask.*[on1-open1]/mapdel/gamma_bar*2*pi;
%Bz_1= [on1-open1]/mapdel/gamma_bar*2*pi;
write_img('Bz1.img',Bz_1,hdr1);
write_hdr('Bz1.hdr',hdr1);

ortho2005([],'anat_file','Bz1.img', 'wscale',[])


