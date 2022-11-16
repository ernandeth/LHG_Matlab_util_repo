
%dbstop if error
close all
clear all

if 1
    totalfig = figure;
    plotfig = figure;
    wirefig = figure;
    Bshieldfig = figure;
    Topshieldfig = figure;
    coilfig = figure;
end

Nvox=64;

zview = Nvox/2 +7; % corresponds to 0.0625 : halfway between coil and head center
yview = Nvox/2;
xview = Nvox/2;

make_masks

if 0
    figure(wirefig)
    ov([],10*(headmask + cylmask +targetmask),...
        xview, yview, zview,0);
end

current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

config='hemisphere'
r = 0.06;  % this radius is a fraction of the FOV
aperture = pi/6;

R = 0.10;
shield_offset = .0625;
theta = pi/8;

FOV = 0.5; % m
% this means that each voxel's size  is 0.5/64 meters = 0.0078 m

% Here are figures of merit from the flat figure 8 coil for comparison in
% the file flat8_measures.mat.
% you can compare these to those.
if 0
    load Efield_fig8_01
    targetN = sum(targetmask(:));
    Ehead = E1 .* headmask;
    flat8_sharpness = sum(E1(:).*cylmask(:)) / sum(Ehead(:))
    flat8_precision =   sum(E1(:).*targetmask(:)) / sum(Ehead(:))
    flat8_targetfield = sum(E1(:).*targetmask(:)) /targetN

    save flat8_measures.mat flat8_precision flat8_sharpness  flat8_targetfield

else
    load flat8_measures.mat
end


if 1

    % Calculate E field from coil alone:
    figure(coilfig)
    hold on
    %[E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, theta);
    %[E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, 0);
    %save Efield_fig8_01 E1 Ex1 Ey1 Ez1
    load Efield_fig8_01.mat

    % big lower shield is #2
    %[E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.07, 0.125 - shield_offset, 0, pi/6, 0);
    %[E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.07, 0.125 - shield_offset, 0, pi/6, 0);
    %save Efield_bottomshield_01 E2 Ex2 Ey2 Ez2
    load Efield_bottomshield_04

    % little top shield is #3
    %[E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/8, 0);
    % [E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/8, 0);
    % save Efield_topshield_01 E3 Ex3 Ey3 Ez3
    load Efield_topshield_03

else
    load Ea3.mat
end

figure(Bshieldfig);
ov([],headmask.*E2 *256/(max(E2(:)) - min(E2(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', 'shield coil')

figure(Topshieldfig);
ov([],headmask.*E3 *256/(max(E3(:)) - min(E3(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', 'shield coil')

figure(coilfig)
ov([],headmask.*E1 *256/(max(E1(:)) - min(E1(:))),...
    xview, yview, zview,0);
set(gcf, 'Name', 'main coil')

drawnow


counter=1;
fwhm = zeros(21);
sharp = zeros(21);
prec = zeros(21);
targetfield = (21);


mysfactorsBottom = linspace(-2, 2,21);   %  lower shield (big)
mysfactorsTop = linspace(-2,2,21);  % upper shield

for sfactorBottom = mysfactorsBottom
    for sfactorTop = mysfactorsTop

        Ex = sfactorBottom*Ex2 + sfactorTop*Ex3 + Ex1;
        Ey = sfactorBottom*Ey2 + sfactorTop*Ey3 + Ey1;
        Ez = sfactorBottom*Ez2 + sfactorTop*Ez3 + Ez1;

        E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);

        if 0
            figure(totalfig)
            ov([],E *256/(max(E(:)) - min(E(:))) ,...
                xview, yview, zview,0);
            set(gcf, 'Name', sprintf('combination %f , %f shielding', sfactorBottom, sfactorTop))
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

        end
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
sharp = reshape(sharp,21,21) / flat8_sharpness;
prec = reshape(prec,21,21)/ flat8_precision;
targetfield = reshape(targetfield, 21,21)/flat8_targetfield;


figure; imagesc(targetfield); axis xy ; title('Field on Target')
ylabel('top shield effect')
set(gca,'ytick',linspace(1,21,5));
set(gca,'yticklabel', linspace(min(mysfactorsTop),max(mysfactorsTop),5))
xlabel('bottom shield effect'),
set(gca,'xtick',linspace(1,21,5));
set(gca,'xticklabel', linspace(min(mysfactorsBottom),max(mysfactorsBottom),5))
colorbar

figure; imagesc(prec); axis xy ; title('precision')
ylabel('top shield effect')
set(gca,'ytick',linspace(1,21,5));
set(gca,'yticklabel', linspace(min(mysfactorsTop),max(mysfactorsTop),5))
xlabel('bottom shield effect'),
set(gca,'xtick',linspace(1,21,5));
set(gca,'xticklabel', linspace(min(mysfactorsBottom),max(mysfactorsBottom),5))
colorbar

figure; imagesc(sharp); axis xy; title('sharpness')
ylabel('top shield effect')
set(gca,'ytick',linspace(1,21,5));
set(gca,'yticklabel', linspace(min(mysfactorsTop),max(mysfactorsTop),5))
xlabel('bottom shield effect'),
set(gca,'xtick',linspace(1,21,5));
set(gca,'xticklabel', linspace(min(mysfactorsBottom),max(mysfactorsBottom),5))
colorbar


% show the best combination
ind = find(sharp(:)==max(sharp(:)));
[i j] = ind2sub(size(sharp),ind);

sfactorTop = mysfactorsTop(i);
sfactorBottom = mysfactorsBottom(j);

Ex = sfactorBottom*Ex2 + sfactorTop*Ex3 + Ex1;
Ey = sfactorBottom*Ey2 + sfactorTop*Ey3 + Ey1;
Ez = sfactorBottom*Ez2 + sfactorTop*Ez3 + Ez1;

E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);


figure(totalfig)
ov([],headmask.*E *256/(max(E(:)) - min(E(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', sprintf('Best combination %f , %f shielding', sfactorBottom, sfactorTop))
p = get(gcf,'Position');
p(1) = p(1)+500;

figure(plotfig);
set(gcf,'Position',p);
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

fprintf('\nMax. sharpness and precision relative to flat fig-8 alone: %f %f \n', max(sharp(:)), max(prec(:)))

save Ea3.mat

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
