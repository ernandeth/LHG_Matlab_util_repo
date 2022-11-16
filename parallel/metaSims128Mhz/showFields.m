
load coilmeta1.mat
Bx1 = abs(Bx); 
By1 = abs(By); 
B1 = complex(Bx1 +By1);

load coilnometa1.mat
Bx1no = abs(Bx); 
By1no = abs(By); 
B1no = complex(Bx1no +By1no);

subplot(211)
imagesc(abs(B1))
subplot(212)
imagesc(abs(B1no))
pause(0.2)


load coilmeta2.mat
Bx2 = abs(Bx); 
By2 = abs(By); 
B2 = complex(Bx2 +By2);

load coilnometa2.mat
Bx2no = abs(Bx); 
By2no = abs(By); 
B2no = complex(Bx2no +By2no);

subplot(211)
imagesc(abs(B2))
subplot(212)
imagesc(abs(B2no))
pause(0.2)



load coilmeta3.mat
Bx3 = abs(Bx); 
By3 = abs(By); 
B3 = complex(Bx3 +By3);

load coilnometa3.mat
Bx3no = abs(Bx); 
By3no = abs(By); 
B3no = complex(Bx3no +By3no);

subplot(211)
imagesc(abs(B3))
subplot(212)
imagesc(abs(B3no))
pause(0.2)




load coilmeta4.mat
Bx4 = abs(Bx); 
By4 = abs(By); 
B4 = complex(Bx4 +By4);

load coilnometa4.mat
Bx4no = abs(Bx); 
By4no = abs(By); 
B4no = complex(Bx4no +By4no);

subplot(211)
imagesc(abs(B4))
subplot(212)
imagesc(abs(B4no))
pause(0.2)




load coilmeta5.mat
Bx5 = abs(Bx); 
By5 = abs(By); 
B5 = complex(Bx5 +By5);

load coilnometa5.mat
Bx5no = abs(Bx); 
By5no = abs(By); 
B5no = complex(Bx5no +By5no);

subplot(211)
imagesc(abs(B5))
subplot(212)
imagesc(abs(B5no))
pause(0.2)




load coilmeta6.mat
Bx6 = abs(Bx); 
By6 = abs(By); 
B6 = complex(Bx6 +By6);

load coilnometa6.mat
Bx6no = abs(Bx); 
By6no = abs(By); 
B6no = complex(Bx6no +By6no);

subplot(211)
imagesc(abs(B6))
subplot(212)
imagesc(abs(B6no))
pause(0.2)


load coilmeta7.mat
Bx7 = abs(Bx); 
By7 = abs(By); 
B7 = complex(Bx7 +By7);

load coilnometa7.mat
Bx7no = abs(Bx); 
By7no = abs(By); 
B7no = complex(Bx7no +By7no);

subplot(211)
imagesc(abs(B7))
subplot(212)
imagesc(abs(B7no))
pause(0.2)




load coilmeta8.mat
Bx8 = abs(Bx); 
By8 = abs(By); 
B8 = complex(Bx8 +By8);

load coilnometa8.mat
Bx8no = abs(Bx); 
By8no = abs(By); 
B8no = complex(Bx8no +By8no);

subplot(211)
imagesc(abs(B8))
subplot(212)
imagesc(abs(B8no))
pause(0.2)


all_Sno = [
    B1no(:)'; 
    B2no(:)';
    B3no(:)';
    B4no(:)';
    B5no(:)';
    B6no(:)';
    B7no(:)';
    B8no(:)';
    ];

all_S = [
    B1(:)'; 
    B2(:)';
    B3(:)';
    B4(:)';
    B5(:)';
    B6(:)';
    B7(:)';
    B8(:)';
    ];

