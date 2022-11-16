%%%
deg=30 
current =1;
config='fig8'

switch (config)
    case 'fig8'
        %deg=30
        ang = deg*pi/180;
        
        wx = 0.6 + 0.3*cos(ang)*sin([0:0.5:2*pi]);
        wy = 0 + 0.3*cos(ang)*cos([0:0.5:2*pi]);
        wz = 0 + 0.3*sin(ang)*sin([0:0.5:2*pi]);
        wire1=[wx' wy' wz'];
        
        wx = -0.6 - 0.3*sin([0.1:0.5:2*pi]);
        wy = 0 - 0.3*cos(ang)*cos([0.1:0.5:2*pi]);
        wz = 0 - 0.3*cos(ang)*sin(ang)*sin([0:0.5:2*pi]);
        
        wire2=[wx' wy' wz'];
        
        [B1 Bx1, By1, Bz1]=biot3d(current, wire1);
        [B2 Bx2, By2, Bz2]=biot3d(current, wire2);
        Bx = Bx1 + Bx2;
        By = By1 + By2;
        Bz = Bz1 + Bz2;
        B=B1+B2;
    case 'helmholtz'
        deg=0;
        ang = 0;
        wx = 0 + 0.3*cos(ang)*sin([0:0.5:2*pi]);
        wy = 0 + 0.3*cos(ang)*cos([0:0.5:2*pi]);
        wz = -0.3 + 0.3*sin(ang)*sin([0:0.5:2*pi]);
        wire1=[wx' wy' wz'];
        
        wx = 0 + 0.3*cos(ang)*sin([0.1:0.5:2*pi]);
        wy = 0 + 0.3*cos(ang)*cos([0.1:0.5:2*pi]);
        wz = 0.3 - 0.3*sin(ang)*sin([0:0.5:2*pi]);
        wire2=[wx' wy' wz'];
        
        [B1 Bx1, By1, Bz1]=biot3d(current, wire1);
        [B2 Bx2, By2, Bz2]=biot3d(current, wire2);
        Bx = Bx1 + Bx2;
        By = By1 + By2;
        Bz = Bz1 + Bz2;
        B=B1+B2;
        
   
end



DIM=length(Bx);


% now add another couple of coils perpendicular...
for z=1:size(Bx,3)
    Bx(:,:,z) = Bx(:,:,z) + Bx(:,:,z)';
    By(:,:,z) = By(:,:,z) + By(:,:,z)';
    Bz(:,:,z) = Bz(:,:,z) + Bz(:,:,z)';
end

%Bmag =sqrt(B(:,1).^2 + B(:,2).^2 + B(:,3).^2);
Bmag =reshape(Bmag,DIM,DIM,DIM);

Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);

yplane = 0
zplane=0.3

zplane = (zplane + 1)/0.1;
yplane = (yplane + 1)/0.1;

figure
subplot(2,2,1)
imagesc(Bz1(:,:,10))
title('coil 1')
subplot(2,2,2)
imagesc( Bz2(:,:,zplane)), title('coil 2')
subplot(2,2, 3)
imagesc( Bz(:,:,zplane)), title('coil 1+2 (Z plane)')
subplot(2,2,4)
imagesc(squeeze(Bz(:,yplane,:))),title('coil 1+2 (Y plane)')
colorbar

figure
subplot(211)
plot3(wire1(:,1),wire1(:,2),wire1(:,3) )
hold on, 
plot3(wire2(:,1),wire2(:,2),wire2(:,3) )
grid on, axis([-1 1 -1 1 -1 1])

subplot(212), hold on
plot([-1:0.1:1], squeeze(Bz(10,10,:)))
plot([-1:0.1:1], squeeze(Bmag(10,10,:)),'g')
xlabel('Z'), ylabel('Bz, ||B||')
legend('Bz', '||B||')

eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bx.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 By.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bz.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bmag.hdr', DIM, DIM, DIM))
h=read_hdr('Bz.hdr');
write_img('Bx.img', 1e10*Bx, h);
write_img('By.img', 1e10*By, h);
write_img('Bz.img', 1e10*Bz, h);
write_img('Bmag.img', 1e10*Bmag, h);

eval(sprintf('! mkdir %ddeg', deg));
eval(sprintf('! mv *.img %ddeg ; mv *.hdr %ddeg', deg, deg));
