
%dbstop if error
close all

totalfig = figure;
plotfig = figure;
wirefig = figure;
shieldfig = figure;
coilfig = figure;

make_masks
ov([],10*(headmask + cylmask +targetmask),...
    xview, yview, zview,0);

current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

config='hemisphere'
r = 0.06;  % this radius is a fraction of the FOV
aperture = pi/6;

R = 0.10;
sfactor = -0.33;
shield_offset = .0625;


theta = pi/8;
%theta = 0;

FOV = 0.5; % m

% Calculate E field from coil alone:
[E1 E1x E1y E1z] =  make_fig8(current, r, 0.125, -theta);
[Es1 Esx1 Esy1 Esz1] = make_doubleC(current, 0.08, 0.125 - shield_offset, 0, pi/6);
[Es2 Esx2 Esy2 Esz2] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/8);


zview = 40; % corresponds to 0.0625 : halfway between coil and head center
yview = 32;
xview = 32;

figure(shieldfig);
ov([],Es1 *256/(max(Es1(:)) - min(Es1(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', 'shield coil')

figure(shieldfig);
ov([],Es2 *256/(max(Es2(:)) - min(Es2(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', 'shield coil')

figure(coilfig)
ov([],E1 *256/(max(E1(:)) - min(E1(:))),...
    xview, yview, zview,0);
set(gcf, 'Name', 'main coil')


counter=1;
fwhm = zeros(21);
sharp = zeros(21);
prec = zeros(21);
targetfield = (21);
for sfactor1 = linspace(-1, 1,21);
    for sfactor2 = linspace(-1,1,21);
        
        Ex = sfactor1*Esx1 + sfactor2*Esx2 + E1x;
        Ey = sfactor1*Esy1 + sfactor2*Esy2 +E1y;
        Ez = sfactor1*Esz1 + sfactor2*Esz2 +E1z;

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
        hold off

        subplot (3,1,2)
        plot(squeeze(E(xview,:,zview))), title('y axis')
        hold on
        plot(squeeze(E1(xview,:,zview)),'k'), title('y axis')
        hold off

        subplot (3,1,3)
        plot(squeeze(E(xview, yview,:))), title('z axis')
        hold on
        plot(squeeze(E1(xview, yview,:)),'k'), title('z axis')
        hold off


        % calculate measures of goodness
        targetN = sum(targetmask(:));
        Ehead = E .* headmask;
        tmp = my_fwhm(Ehead(:,:,zview)); tmp = fwhm(counter);
        sharp(counter) = sum(E(:).*cylmask(:)) / sum(Ehead(:));
        prec(counter) =  sum(E(:).*targetmask(:)) / sum(Ehead(:));
        targetfield(counter) = sum(E(:).*targetmask(:)) /targetN ;
        counter = counter+1;
    end
end

fwhm = reshape(fwhm, 21,21);
sharp = reshape(sharp,21,21);
prec = reshape(prec,21,21);
targetfield = reshape(targetfield, 21,21);


figure; imagesc(targetfield); axis xy
figure; imagesc(prec); axis xy       
figure; imagesc(sharp); axis xy
xlabel('top shield effect'), ylabel('bottom shield effect')
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
