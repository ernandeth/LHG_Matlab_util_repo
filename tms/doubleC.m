%dbstop if error
%close all
doBatch = 1;

if ~doBatch
% this is the stuff we do one time only)
plotfig = figure;
wirefig = figure;
coilfig = figure;
roifig = figure;

deg=0
current = 1e3/1e-4;  % this is actually current ramp rate.  
               % We get to 1 KAmp in 0.1 ms

config='hemisphere'
shield_offset=0;

R = 0.08; % m
FOV = 0.5;

wire=[];
axis xy

zview = 32; % corresponds to 0.25 : halfway between coil and head center
yview = 32;
xview = 32;

theta = pi/6;
theta = 0;

aperture = pi/6;

end

loopSegs = linspace( 2*aperture, 2*pi-2*aperture, 20);

wx = R*sin(loopSegs);
wy = -R*cos(loopSegs);
wz =  0 -shield_offset + wy*sin(theta);

wire =[
    0 0 10
    wx' R+wy' wz';
	0 0 10;
	wx' -R-wy' wz';
    0 0 10.01
    ];

figure(wirefig)
plot3(wire(:,1), wire(:,2), wire(:,3),'b');
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
drawnow

% Allocate space and define some conductivity space
Nvox = 64;

% Calculate E field from coil alone:
[E Ex Ey Ez] =  Efield(current, wire, [Nvox Nvox Nvox], FOV);
E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);

% make a headsphere
headmask = zeros(size(E));
center = [32,32,30];
head_radius = 16;
for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
			if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= head_radius,
				headmask(x,y,z) = 1;
			end,
		end,end,end


% make a target
targetmask = zeros(size(E));
center = [32,32,30];
target_radius = 2;
for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
			if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= target_radius,
				targetmask(x,y,z) = 1;
			end,
		end,end,end


% make a cylinder
cylmask = zeros(size(E));
head_radius = 16;
center = [32,32,30];
for x=1:Nvox, for y=1:Nvox, for z=center-head_radius-1:center+head_radius-1
			if sqrt((x-center(1))^2 +(y-center(2))^2 ) <= 1,
				cylmask(x,y,z) = 1;
			end,
		end,end,end


figure(roifig)
ov([],10*(headmask + cylmask +targetmask),...
	xview, yview, zview,0);


figure(coilfig)
ov([],E * 256/(max(E(:))-min(E(:))),xview, yview, zview,0);
set(gcf, 'Name', 'main coil')

p = get(gcf,'Position');
p(1) = p(1)+500;

xx = linspace(-FOV/2, FOV/2, 64);
figure(plotfig); set(gcf,'Position',p);
normE = E / sum(E(:));

subplot (3,1,1)
plot(xx,squeeze(E(:,yview, zview))), title('x axis')
fatlines, dofontsize(16);
xlabel('distance (m)')
ylabel('E field (V/m)');

hold on

subplot (3,1,2)
plot(xx,squeeze(E(xview,:,zview))), title('y axis')
fatlines, dofontsize(16);
xlabel('distance (m)')
ylabel('E field (V/m)');
hold on

subplot (3,1,3)
plot(xx,squeeze(E(xview, yview,:))), title('z axis')
fatlines, dofontsize(16);
xlabel('distance (m)')
ylabel('E field (V/m)');
hold on



save C_coil.mat

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zview = 60
figure
subplot(211)
quiver(squeeze(Esx(1:5:end,1:5:end,zview)), squeeze(Esy(1:5:end,1:5:end,zview)),'r');
axis square; axis xy; axis tight; title(sprintf('shield slice %d transverse', zview));

subplot(212)
quiver(squeeze(Esx(1:5:end, yview, 1:3:end)), squeeze(Esz(1:5:end, yview, 1:3:end)));
axis square; axis xy; axis tight; title(sprintf('shield slice %d coronal', yview));

figure
subplot(211)
quiver(squeeze(E1x(1:5:end,1:5:end,zview)), squeeze(E1y(1:5:end,1:5:end,zview)),'r');
axis square; axis xy; axis tight; title(sprintf('fig8 slice %d transverse', zview));

subplot(212)
quiver(squeeze(E1x(1:5:end, yview, 1:3:end)), squeeze(E1z(1:5:end, yview, 1:3:end)));
axis square; axis xy; axis tight; title(sprintf('fig8 slice %d coronal', yview));

figure
subplot(211)
quiver(squeeze(Ex(1:5:end,1:5:end,zview)), squeeze(Ey(1:5:end,1:5:end,zview)));
axis square; axis xy; axis tight; title(sprintf('total slice %d transverse', zview));

subplot(212)
quiver(squeeze(Ex(1:5:end, yview, 1:5:end)), squeeze(Ez(1:5:end, yview, 1:5:end)));
axis square; axis xy; axis tight; title(sprintf('total slice %d coronal', yview));
