%%%
deg=0
current =1;  % Amp
config='hemisphere'
r = 0.15;
wire=[];

angle = theta;

wx = 0 + r*cos(angle)*sin([0:0.5:2*pi]);
wy = 0 +            r*cos([0:0.5:2*pi]);
wz = 0*ones(size(wx));
wire=[wire;  wx' wy' wz'];
wire = [wire; wire(1,:)];

plot3(wire(:,1),wire(:,2),wire(:,3) )
grid on, axis([-1 1 -1 1 -1 1])


% do the biot-savart law here
[B Bx, By, Bz]=biot3d(current, wire);
Bz_old = Bz;
DIM=length(Bx);

% add a 0.6 G/cm gradient to the whole thing
[x y Bzgrad] = meshgrid(zeros(1,DIM), zeros(1,DIM), linspace(-100*0.6,100*0.6, DIM));

Bz = 1e7*Bz_old+ Bzgrad;

Bmag = (Bx.^2 + By.^2 + (Bz).^2).^0.5;

xplane = floor(DIM/2);
yplane = floor(DIM/2);
zplane = floor(DIM/2);

for n=1:DIM
    ov([],(Bz), xplane,yplane,n,0);
    drawnow; pause(0.1)
end



lightbox(-Bz_old, [-1e-5 1e-5],10);

eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bx.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 By.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bz.hdr', DIM, DIM, DIM))
eval(sprintf('! avwcreatehd %d %d %d 1 1 1 1 1 0 0 0 4 Bmag.hdr', DIM, DIM, DIM))
h=read_hdr('Bz.hdr');
write_img('Bx.img', 1e10*Bx, h);
write_img('By.img', 1e10*By, h);
write_img('Bz.img', 1e10*Bz, h);
write_img('Bmag.img', 1e10*Bmag, h);
