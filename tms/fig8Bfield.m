% Simplefigure 8 magnetic field
wire =[];
r = 0.05;
theta=70*2*pi/180;
zpos=0;
xgap=0;

% Make the fig8 coils:
loopSegs = linspace(0,2*pi,15);

wx = r*sin(loopSegs) ;
wy = -r*cos(loopSegs);
wz = zpos - wy*sin(theta) - r*sin(theta) ;

loop =[wx'  -r-wy'-xgap/2 wz';
    %0 0 0;
    wx' r+wy'+xgap/2 wz'];

wire = [wire; loop];

current = 6* 14/50; % 6 turns x 14 V / 50 ohms

% this produces a field with FOV=0.05 m, ad 0.01 cm resolution
[B, Bx, By , Bz]=biot3d(current, wire)

subplot(211), plot3(wire(:,1), wire(:,2), wire(:,3))

subplot(212),    lightbox(B);

plot([-0.1:0.01:0.1], squeeze(6*B( 5,11,:)) ) ;


