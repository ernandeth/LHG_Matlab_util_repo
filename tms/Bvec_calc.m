cd open1
! gsp21a -m P*
! ref-stack open1 64 64 30 1 1 1 ref.64.s
loc1_open = read_img2('open1.img');
cd ..

cd on1
! gsp21a -m P*
! ref-stack on1 64 64 30 1 1 1 ref.64.s
loc1_on = read_img2('on1.img');
figure
lightbox(loc1_on-loc1_open); title('Phase change  position 1')
cd ..


cd open2
! gsp21a -m P*
! ref-stack open2 64 64 30 1 1 1 ref.64.s
loc2_open = read_img2('open2.img');
cd ..



cd on2
! gsp21a -m P*
! ref-stack on2 64 64 30 1 1 1 ref.64.s
loc2_on = read_img2('on2.img');
figure
lightbox(loc2_on-loc2_open); title('Phase change  position 1')
cd ..


cd open3
! gsp21a -m P*
! ref-stack open3 64 64 30 1 1 1 ref.64.s
loc3_open = read_img2('open3.img');
cd ..



cd on3
! gsp21a -m P*
! ref-stack on3 64 64 30 1 1 1 ref.64.s
loc3_on = read_img2('on3.img');
figure
lightbox(loc3_on-loc3_open); title('Phase change  position 1')
cd ..

!ln -s */.img .
!ln -s */.hdr .
!avwmerge -t tmp open*.img

