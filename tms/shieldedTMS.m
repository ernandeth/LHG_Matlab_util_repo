
%dbstop if error
close all
doBatch = 0

if ~doBatch  % simulate the TMS coil alone
    % this is the stuff we do one time only)
    totalfig = figure;
    plotfig = figure;
    wirefig = figure;
    shieldfig = figure;
    coilfig = figure;

    deg=0
    current =1e3/1e-4;
    config='hemisphere'
    r = 0.06;  % this radius is a fraction of the FOV


    R = 0.10;
    sfactor = -0.33;
    shield_offset = .0625;

    
    wire=[];
    axis xy

    zview = 40; % corresponds to 0.0625 : halfway between coil and head center
    yview = 32;
    xview = 32;

    theta = pi/8;
    theta = 0;

    FOV = 0.5; % m

    % Make the fig8 coils:
    loopSegs = linspace(0,2*pi,20);
    %r = 0.2;
    wx = r*sin(loopSegs);
    wy = -r*cos(loopSegs);
    wz = 0.125 + wy*sin(theta);
    % in the flat case, I want the apex of the coil in the same location as the
    % bent one
    %wz(:) = 0.5 + max(wy)*sin(pi/6);

    %loop =[wx' wy' wz';];
    loop =[wx' -r-wy' wz';
        wx' r+wy' wz'];

    wire = [wire; loop];

    figure(wirefig)
    plot3(wire(:,1), wire(:,2), wire(:,3),'k')
    %hold on,plot(wire(:,1), wire(:,2),'r')
    hold on , plot3(0,0,0,'ro')
    axis([-1 1 -1 1 -1 1])
    %plot3(wire(1,1), wire(1,2), wire(1,3),'g*')
    %plot3(wire(round(end/2),1), wire(round(end/2),2), wire(round(end/2),3),'go')
    %plot3(wire(end,1), wire(end,2), wire(end,3),'g+')
    drawnow
    % Allocate space and define some conductivity space
    Nvox = 64;

    % Calculate E field from coil alone:
    [E1 E1x E1y E1z] =  Efield(current, wire, floor([Nvox Nvox Nvox]),FOV);

    % make a headsphere
    headmask = zeros(size(E1));
    center = [32,32,32];
    head_radius = 16;
    for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
                if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= head_radius,
                    headmask(x,y,z) = 1;
                end,
            end,end,end


    % make a target
    targetmask = zeros(size(E1));
    center = [32,32,40];
    target_radius = 2;
    for x=1:Nvox, for y=1:Nvox, for z=1:Nvox
                if sqrt((x-center(1))^2 +(y-center(2))^2 + (z-center(3))^2) <= target_radius,
                    targetmask(x,y,z) = 1;
                end,
            end,end,end


    % make a cylinder
    cylmask = zeros(size(E1));
    head_radius = 16;
    center = [32,32,30];
    for x=1:Nvox, for y=1:Nvox, for z=center-head_radius-1:center+head_radius-1
                if sqrt((x-center(1))^2 +(y-center(2))^2 ) <= 1,
                    cylmask(x,y,z) = 1;
                end,
            end,end,end


    figure()
    ov([],10*(headmask + cylmask +targetmask),...
        xview, yview, zview,0);


end


figure(wirefig);
plot3(wire(:,1), wire(:,2), wire(:,3),'k'); hold on

% Now simulate the shield:
aperture = pi/6;
loopSegs = linspace( aperture, 2*pi - aperture, 20);

wx = R*sin(loopSegs);
wy = -R*cos(loopSegs);
wz = 0.125 - shield_offset + wy*sin(theta);

swire =[0 0 10;
    wx' R+wy' wz';
    0 0 10;
    wx' -R-wy' wz';
    0 0 10.01];

plot3(swire(:,1), swire(:,2), swire(:,3),'b');
axis([-1 1 -1 1 -1 1]/2);

drawnow
%
% % for testing purposes
% wire =[wx' r+wy' wz'];
% plot3(wire(1:end-5,1), wire(1:end-5,2), wire(1:end-5,3),'k')
% hold on %, plot3(0,0,0,'ro')
% wire =[wx' -r-wy' wz'];
% plot3(wire(:,1), wire(:,2), wire(:,3),'k')

% plot3(wire(1,1), wire(1,2), wire(1,3),'g*')
% plot3(wire(round(end/2),1), wire(round(end/2),2), wire(round(end/2),3),'go')
% plot3(wire(end,1), wire(end,2), wire(end,3),'g+')



[Es Esx Esy Esz] =  Efield(current, swire, [Nvox Nvox Nvox], FOV);

figure(shieldfig);
ov([],Es *256/(max(Es(:)) - min(Es(:))) + 0*(headmask +cylmask +targetmask),...
    xview, yview, zview,0);
set(gcf, 'Name', 'shield coil')

figure(coilfig)
ov([],E1 *256/(max(E1(:)) - min(E1(:))) + 0*(headmask +cylmask +targetmask),...
    xview, yview, zview,0);
set(gcf, 'Name', 'main coil')

Ex = sfactor*Esx + E1x;
Ey = sfactor*Esy + E1y;
Ez = sfactor*Esz + E1z;

E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);

figure(totalfig)
ov([],E *256/(max(E(:)) - min(E(:))) + 0*(headmask +cylmask +targetmask),...
    xview, yview, zview,0);
set(gcf, 'Name', sprintf('combination %f shielding', sfactor))
p = get(gcf,'Position');
p(1) = p(1)+500;

figure(plotfig); set(gcf,'Position',p);
normE = E / sum(E(:));
subplot (3,1,1)
plot(squeeze(E(:,yview, zview))), title('x axis')
hold on
plot(squeeze(E1(:,yview, zview)),'k'), title('x axis')

subplot (3,1,2)
plot(squeeze(E(xview,:,zview))), title('y axis')
hold on
plot(squeeze(E1(xview,:,zview)),'k'), title('y axis')

subplot (3,1,3)
plot(squeeze(E(xview, yview,:))), title('z axis')
hold on
plot(squeeze(E1(xview, yview,:)),'k'), title('z axis')


% calculate measures of goodness
targetN = sum(targetmask(:));
Ehead = E .* headmask;
fwhm = my_fwhm(Ehead(:,:,zview))
sharp = sum(E(:).*cylmask(:)) / sum(Ehead(:))
prec =  sum(E(:).*targetmask(:)) / sum(Ehead(:))
targetfield = sum(E(:).*targetmask(:)) /targetN


save Ea.mat

return

figure(totalfig)
Ediff=E1-E;
ov([],Ediff *8*512/(max(Ediff(:)) - min(Ediff(:))) + 0*(headmask +cylmask +targetmask),...
    xview, yview, zview,0);
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
