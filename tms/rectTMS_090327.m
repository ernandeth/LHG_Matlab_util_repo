%dbstop if error
close all
deg=0 
current =1;
config='rect'
omega = 10000; %hz

wire = [0 0.05  0 ;
        0 -0.05 0 ;
        0 -0.05 0.5;
        0 0.05  0.5;
        0 0.05 0];
    

plot3(wire(:,1), wire(:,2), wire(:,3),'k')
hold on , plot3(0,0,0,'ro')
axis([-1 1 -1 1 -1 1])

wire2 = [0.05 0.05  0 ;
        0.05 -0.05 0 ;
        0.08 -0.05 0.5;
        0.08 0.05  0.5;
        0.05 0.05 0];

plot3(wire2(:,1), wire2(:,2), wire2(:,3),'k')

wire3 = [-0.05 0.05  0 ;
        -0.05 -0.05 0 ;
        -0.08 -0.05 0.5;
        -0.08 0.05  0.5;
        -0.05 0.05 0];

plot3(wire3(:,1), wire3(:,2), wire3(:,3),'k')

drawnow
% Allocate space and define some conductivity space
Nvox = 64;
zview = 30;

% Calculate E field from coil alone:    


[E1 E1x E1y E1z] =  Helmholtz3d(current, omega, wire);

[E2 E2x E2y E2z] =  Helmholtz3d(current, omega, wire2);

[E3 E3x E3y E3z] =  Helmholtz3d(current, omega, wire3);

Ex = imag(E1x + E2x + E3x);
Ey = imag(E1y + E2y + E3y);
Ez = imag(E1z + E2z + E3z);

Enorm = sqrt(Ex.^2 +Ey.^2 + Ez.^2);

figure

for c=15:35
    ov([],(Ey),21 ,21,c,0);
    pause(0.2)
end
   

%save Ea.mat 

return

zview = 20
yview = 20;

figure
subplot(221)

quiver(squeeze(Ex(:,:,zview)), squeeze(Ey(:,:,zview)));

axis xy; title(sprintf('slice %d transverse XY', zview));

subplot(222)
quiver(squeeze(Ex(1:5:end, yview, 1:3:end)), squeeze(E1z(1:5:end, yview, 1:3:end)));
axis xy; title(sprintf('slice %d coronal XZ', yview));

figure
subplot(223)
quiver(squeeze(Ex(1:5:end,1:5:end,zview)), squeeze(Ey(1:5:end,1:5:end,zview)));
axis xy; title(sprintf('total slice %d transverse', zview));

subplot(212)
quiver(squeeze(Ex(1:5:end, yview, 1:5:end)), squeeze(Ez(1:5:end, yview, 1:5:end)));
axis xy; title(sprintf('total slice %d coronal', yview));
