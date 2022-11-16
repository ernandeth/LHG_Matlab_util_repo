function [cost, E] = shieldedTMS_lsq(parms, consts)
%
% function [cost, E] = shieldedTMS_lsq(parms, consts)
%
% sfactorBottom = parms(1);
% sfactorTop = parms(2);
% theta = parms(3);
% R2 = parms(4);
% R3 = parms(5);
% Zbottom = parms(6);
% top_Z = parms(7);
%
% consts = [] for now
%
% Note: if the doDivergence=1, then we optimize on the divergence
% of the E field rather than on the electric field
parms

doDivergence=1;

% sfactorBottom = parms(1);
% sfactorTop = parms(2);
% theta = parms(3);
% Rtop = parms(4);
% Rbottom = parms(5);
% Zbottom = parms(6);
% top_Z = parms(7);

sfactorBottom = parms(1);
xgap = parms(2);
theta = parms(3);
Rbottom = parms(4);
Zbottom = parms(5);


FOV = 0.5; % m
Nvox=33;
zview = round(Nvox/2 + Nvox/8); % corresponds to 0.0625 : halfway between coil and head center
yview = round(Nvox/2);
xview = round(Nvox/2);

make_masks


current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

R1 = 0.06;  % m
aperture = pi/6;



% Calculate E field from coil alone:
[E1 Ex1 Ey1 Ez1] =  make_fig8(current, R1, xgap, 0.125, theta);
% bottom shield
[E2 Ex2 Ey2 Ez2] = make_doubleC(current, Rbottom, Zbottom, 0, pi/6, 0);
% top shield
%[E3 Ex3 Ey3 Ez3] = make_doubleC(current, Rtop, top_Z, 0, pi/8, 0);

% Ex = sfactorBottom*Ex2 + sfactorTop*Ex3 + Ex1;
% Ey = sfactorBottom*Ey2 + sfactorTop*Ey3 + Ey1;
% Ez = sfactorBottom*Ez2 + sfactorTop*Ez3 + Ez1;
Ex = sfactorBottom*Ex2  + Ex1;
Ey = sfactorBottom*Ey2  + Ey1;
Ez = sfactorBottom*Ez2  + Ez1;

if doDivergence
    mygrid = linspace(-FOV/2,FOV/2, Nvox);
    [x,y,z] = meshgrid(mygrid,mygrid, mygrid ) ;
    E = abs(divergence(x,y,z, (Ex), (Ey), (Ez)));
    E1 = abs(divergence(x,y,z, (Ex1), (Ey1), (Ez1)));
else
    E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);
end
% calculate measures of cost
targetN = sum(targetmask(:));

Ehead = E .* headmask;
Ecyl = E .* cylmask;
sharp = sum(Ecyl(:)) / sum(Ehead(:))
prec =  sum(E(:).*targetmask(:)) / sum(Ehead(:));

% sfigure out penetration within cylinder of interest:
% i.e.:  where is abs(E) = 0.5*abs(E_surface) ?
center =ceil([Nvox, Nvox, Nvox]/2);

Ecyl_max = max(Ecyl(:));
active = zeros(size(Ecyl));
active(Ecyl > Ecyl_max/2) = 1;
inds = find(active);
[x, y, z ]= ind2sub([Nvox, Nvox, Nvox], inds);
z50 = min(z);
z100 = max(z);
[z50 z100]

Esurf = E(center(1), center(2), z100);
E1surf = E1(center(1), center(2), z100);

penetration = E(center(1),center(2), center(3)) / Esurf

%cost = norm([1/sharp  5/penetration])
cost = 1/sharp;

hold off
E(E==max(E(:)))= 0;
E1(E1==max(E1(:)))= 0;


%
% E(E < 0.5*E(center(1), center(2), z100)) =0;
% E1(E1 < 0.5*E1(center(1), center(2), z100)) =0;
%

E = E/Esurf;
E1 = E1/E1surf;
if 1

	sfigure(1);
	ov([], (1+headmask+cylmask)*5 + E *128/(max(E(:)) - min(E(:))) ,...
	    xview, yview, zview,0);
	set(gcf, 'Name', sprintf('shield= %f pen= %f, sharp= %f', sfactorBottom, penetration, sharp))
	drawnow



    p = get(gcf,'Position');
    p(1) = p(1)+500;

    sfigure(2); set(gcf,'Position',p);
    normE = E / sum(E(:));
    subplot (3,1,1)
    plot(squeeze(E(:,yview, zview))), title('x axis')
    hold on
    plot(squeeze(E1(:,yview, zview)),'r'), title('x axis')
    hold off

    subplot (3,1,2)
    plot(squeeze(E(xview,:,zview))), title('y axis')
    hold on
    plot(squeeze(E1(xview,:,zview)),'r'), title('y axis')
    hold off

    subplot (3,1,3)
    plot(squeeze(E(xview, yview,:))), title('z axis')
    hold on
    plot(squeeze(E1(xview, yview,:)),'r'), title('z axis')
    hold off
	sfigure(1);
end
return


