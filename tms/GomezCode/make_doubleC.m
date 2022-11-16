function [wire Jvec] = make_doubleC(current, R, zpos, theta, aperture)
%function [Es Esx Esy Esz] = doubleC(current, R, zpos, bend_theta, aperture)

FOV = 0.5;
wire=[];
Nvox = 64;

loopSegs = linspace( 2*aperture, 2*pi-2*aperture, 20);

wx = R*sin(loopSegs);
wy = -R*cos(loopSegs);
wz =  zpos + wy*sin(theta);

wire =[
    0 0 10
    wx' R+wy' wz';
    0 0 10];
wire=[wire; wx' -R-wy' wz';
    0 0 10.01
    ];
l=length(wire(:,1));
Jvec(l,3)=0;
for i=1:l
    Jvec(i,:)=wire(1+i*(i~=l),:)-wire(i,:);
end



return







