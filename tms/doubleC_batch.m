load C_coil.mat


FWHM=[];
sharpness = [];
precision = [];
field = [];

for aperture=linspace(pi/12,pi/4,5)
    aperture
    for R=linspace(0.06, 0.15,5)
        doubleC
        Ehead = E .* headmask;

        % calculate measures of goodness
        f = my_fwhm(Ehead(:,:,zview));
        FWHM = [FWHM f(1) ];
        sharpness = [sharpness  sum(E(:).*cylmask(:)) / sum(Ehead(:))];
        precision = [precision  sum(E(:).*targetmask(:)) / sum(Ehead(:))];
        field = [field  sum(E(:).*targetmask(:)) ];


        % make some figures
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

        % save things just in case
        save mystate_C.mat  R aperture FWHM sharpness precision field
        %end
    end
end

maskN = sum(targetmask(:));

% new dimensions:  Nsfactor x Nshield_offset x Nradius
FWHM = reshape(FWHM,5,5);
sharpness = reshape(sharpness,5,5);
precision = reshape(precision,5,5);
field = reshape(field, 5,5) / maskN;


aperture=linspace(pi/12,pi/4,5)
R=linspace(0.06, 0.15,5)

tmp = zeros(5,5);


subplot(2,2,1)
lightbox(FWHM), title('FWHM')
ylabel('shield factor')
xlabel('shield radius');
subplot(2,2,2),
lightbox(sharpness), title('sharpness')
subplot(2,2,3)
lightbox(precision), title('precision')
subplot(2,2,4)
lightbox(field), title('Mean Field at Target')


figure
subplot(221)
plot(R, FWHM(:,3)), title('FWHM')
xlabel('Radius')
ylabel('FWHM')

subplot(222)
plot(R, sharpness(:,3)), title('sharpness')
xlabel('Radius')
ylabel('Sharpness')

subplot(223)
plot(R, precision(:,3)), title('precision')
xlabel('Radius')
ylabel('Precision')

subplot(224)
plot(R, field(:,3)), title('E Field')
xlabel('Radius')
ylabel('Mean Field at Target')

figure
subplot(221)
plot(aperture, FWHM(3,:)), title('FWHM')
xlabel('aperture')
ylabel('FWHM')

subplot(222)
plot(aperture, sharpness(3,:)), title('sharpness')
xlabel('aperture')
ylabel('Sharpness')

subplot(223)
plot(aperture, precision(3,:)), title('precision')
xlabel('aperture')
ylabel('Precision')

subplot(224)
plot(aperture, field(3,:)), title('Mean Field at Target')
xlabel('aperture')
ylabel('Field over Target')
