function [Es Esx Esy Esz] = make_doubleC(current, R, zpos, theta, aperture, doSwapxy) 
%function [Es Esx Esy Esz] = doubleC(current, R, zpos, bend_theta, aperture) 
    

FOV = 0.5;
wire=[];
Nvox = 33;

loopSegs = linspace( 2*aperture, 2*pi-2*aperture, 20);

wx = R*sin(loopSegs);
wy = -R*cos(loopSegs);
wz =  zpos + wy*sin(theta);

wire =[
    0 0 10
    wx' R+wy' wz';
	0 0 10;
	wx' -R-wy' wz';
    0 0 10.01
    ];

%%%
if doSwapxy
    tmp=wire;
    tmp(:,1) = wire(:,2);
    tmp(:,2) = wire(:,1);
    wire=tmp;
end
%%
subplot(2,2,4)
plot3(wire(:,1), wire(:,2), wire(:,3),'b');
axis([-0.25 0.25 -0.25 0.25 -0 0.25]);
hold on
drawnow


% Calculate E field from coil alone:
[Es Esx Esy Esz] =  Efield(current, wire, [Nvox Nvox Nvox], FOV);



return



