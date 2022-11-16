%%%
deg=0 
current =1;
config='hemisphere'
r = 0.09
R = 0.6
wire=[];

for theta= 0 : pi/8 : pi-pi/8
    angle = theta;
    
    wx = R*sin(theta) + r*cos(angle)*sin([0:0.5:2*pi]);
    wy = 0 +            r*cos([0:0.5:2*pi]);
    wz = -R*cos(theta)+ r*sin(angle)*sin([0:0.5:2*pi]);
    wire=[wire;  wx' wy' wz']; 
    
end
for theta= pi : pi/8 : 2*pi-pi/8
    angle = theta;
    
    wx = R*sin(theta) - r*cos(angle)*sin([0:0.5:2*pi]);
    wy = 0            - r*cos([0:0.5:2*pi]);
    wz = -R*cos(theta) - r*sin(angle)*sin([0:0.5:2*pi]);
    wire=[wire;  wx' wy' wz']; 
    
end

% add more major loops by rotating the wire over the Z-axis
% for rot=pi/4:pi/4:pi/4
%     mtx=[ 
%         cos(rot)  -sin(rot) 0 ;
%         sin(rot)   cos(rot) 0;
%         0 0 1;
%         ];
%     wire2=mtx*wire';
%     wire = [wire ; wire2'];
% end

% for count=1:length(wire)
%     plot3(wire(1:count,1), wire(1:count,2), wire(1:count,3))
%     axis([-1 1 -1 1 -1 1])
%     drawnow
% end

fprintf('\n Press enter to continue ')
plot3(wire(:,1),wire(:,2),wire(:,3) )
grid on, axis([-1 1 -1 1 -1 1])
pause

% do the biot-savard law here
[B Bx, By, Bz]=biot3d(current, wire);

DIM=length(Bx);


% % % now add another set of coils perpendicular...
% for z=1:size(Bx,3)
%     Bx(:,:,z) = Bx(:,:,z) + By(:,:,z)';
%     By(:,:,z) = By(:,:,z) + Bx(:,:,z)';
%     Bz(:,:,z) = Bz(:,:,z) + Bz(:,:,z)';
% end

Bmag =sqrt(B(:,1).^2 + B(:,2).^2 + B(:,3).^2);
Bmag =reshape(Bmag,DIM,DIM,DIM);

Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);

yplane = 0
zplane=0.3
xplane=0;

zplane = (zplane + 1)/0.1;
xplane = (xplane + 1)/0.1;
yplane = (yplane + 1)/0.1;

figure
subplot(2,2,1)
imagesc( Bmag(:,:,zplane)), title('|B| (Z plane)')
colorbar
subplot(2,2,2)
imagesc(squeeze(Bmag(xplane,:,:))),title('|B|(X plane)')
colorbar
subplot(2,2,3)
imagesc(squeeze(Bmag(:,yplane,:))),title('|B|(y plane)')
colorbar


figure
subplot(2,1,1)
contour( Bmag(:,:,zplane)), title(' (Z plane)')
subplot(2,1,2)
contour(squeeze(Bmag(xplane,:,:)),20),title('(X plane)')
hold on
plot(11+ 10*wire(:,3),11+10*wire(:,1),'k')

figure
plot3(wire(:,1),wire(:,2),wire(:,3) )
grid on, axis([-1 1 -1 1 -1 1])
xlabel('X'), ylabel('Y'), zlabel('Z')

figure
subplot(311)
plot([-1:0.1:1], squeeze(Bmag(:,9:12,11)'))
axis([-1 1 0 1.5e-6])
xlabel('x')

subplot(312)
plot([-1:0.1:1], squeeze(Bmag(9:12,:, 11))')
axis([-1 1 0 1.5e-6])
xlabel('y')

subplot(313), hold on
plot([-1:0.1:1], squeeze(Bmag(11,9:12,:)))
axis([-1 1 0 1.5e-6])
xlabel('Z'), ylabel('||B||')

eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bx.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 By.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bz.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bmag.hdr', DIM, DIM, DIM))
h=read_hdr('Bz.hdr');
write_img('Bx.img', 1e8*Bx, h);
write_img('By.img', 1e8*By, h);
write_img('Bz.img', 1e8*Bz, h);
write_img('Bmag.img', 1e8*Bmag, h);

