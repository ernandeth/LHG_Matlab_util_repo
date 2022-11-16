% call shieldedTMS results and try out combinations of
% shield powers (sfactor)
load Ea.mat
close all


%% iterate over shield sizes
FWHM=[];
sharpness = [];
precision = [];
field=[];

doBatch=1;
for R = linspace(0.08,0.25,10)    
    for shield_offset=linspace(0,0.125,3)
        shieldedTMS
        close(plotfig)
        for sfactor=linspace(-1,1,10)


            % total field:
            Ex = sfactor*Esx + E1x;
            Ey = sfactor*Esy + E1y;
            Ez = sfactor*Esz + E1z;

            E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);
            Ehead = E .* headmask;

            % calculate measures of goodness
            f = my_fwhm(Ehead(:,:,zview));
            FWHM = [FWHM f(1) ];
            sharpness = [sharpness  sum(E(:).*cylmask(:)) / sum(Ehead(:))];
            precision = [precision  sum(E(:).*targetmask(:)) / sum(Ehead(:))];
            field = [field  sum(E(:).*targetmask(:)) ];

            % make some figures
            figure(shieldfig);
            ov([],Es *256/(max(Es(:)) - min(Es(:))) + 0*(headmask +cylmask +targetmask),...
                xview, yview, zview,0);
            set(gcf, 'Name', 'shield coil')

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
            fatlines, dofontsize(16);
            hold on
            subplot (3,1,2)
            plot(squeeze(E(xview,:,zview))), title('y axis')
            fatlines, dofontsize(16);
            hold on
            subplot (3,1,3)
            plot(squeeze(E(xview, yview,:))), title('z axis')
            fatlines, dofontsize(16);
            hold on
            drawnow
        end
        % save things just in case
        save mystateShield.mat field sfactor R shield_offset FWHM sharpness precision
    end
end

%FWHM = FWHM(1:2:end);

% new dimensions:  Nsfactor x Nshield_offset x Nradius
FWHM = reshape(FWHM,10,3,10);
sharpness = reshape(sharpness,10,3,10);
precision = reshape(precision,10,3,10);
field = reshape(field,10,3,10);

sfactors = linspace(-1,1,10)
shield_offsets=linspace(0,0.5,3)
Rs = linspace(0.08,0.25,10)   ;

tmp = zeros(10,10,3);

% new dimensions:  Nsfactor  x Nradius x Nshield_offset
for r=1:3,	tmp(:,:,r) = sharpness(:,r,:);  end
sharpness = tmp;
for r=1:3,	tmp(:,:,r) = FWHM(:,r,:);  end
FWHM = tmp;
for r=1:3,	tmp(:,:,r) = precision(:,r,:);  end
precision = tmp;
for r=1:3,	tmp(:,:,r) = field(:,r,:);  end
field = tmp;


subplot(2,2,1)
lightbox(FWHM), title('FWHM')
ylabel('shield factor')
xlabel('shield radius');
subplot(2,2,2),
lightbox(sharpness), title('sharpness')
subplot(2,2,3)
lightbox(precision), title('precision')
subplot(2,2,4)
lightbox(field), title('Field at Target')


figure
subplot(221)
plot(sfactors, FWHM(:,:,2)), title(['FWHM, z-offset: ' num2str(shield_offsets(3))])
xlabel('Shield Factor')
ylabel('FWHM')

subplot(222)
plot(sfactors, sharpness(:,:,2)), title(['sharpness, z-offset: ' num2str(shield_offsets(3))])
xlabel('Shield Factor')
ylabel('Sharpness')

subplot(223)
plot(sfactors, precision(:,:,2)), title(['Precision, z-offset: ' num2str(shield_offsets(3))])
xlabel('Shield Factor')
ylabel('Precision')

subplot(224)
plot(sfactors, field(:,:,2)), title(['Field at Target, z-offset: ' num2str(shield_offsets(3))])
xlabel('Shield Factor')
ylabel('Field over Target')

