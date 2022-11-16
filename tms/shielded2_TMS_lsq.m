function [cost] = shielded2_TMS_lsq(parms1, parms2, consts)
% 
% function [cost] = shieldedTMS_lsq(parms, consts)
% 
% sfactorBottom = parms(1);
% sfactorTop = parms(2);
% theta = parms(3);
% R2 = parms(4);
% R3 = parms(5);
% Zbottom = parms(6);
% Ztop = parms(7);
%
% consts = [] for now
%


% sfactorBottom = parms(1);
% sfactorTop = parms(2);
% theta = parms(3);
% Rtop = parms(4);
% Rbottom = parms(5);
% Zbottom = parms(6);
% Ztop = parms(7);

sfactorBottom = parms1(1);
xgap = parms1(2);
theta = parms1(3);
Rbottom = parms1(4);
Zbottom = parms1(5);

sfactorTop = parms2(1);
Rtop = parms2(4);
Ztop = parms2(5);

Nvox=33;
zview = round(Nvox/2 + Nvox/8); % corresponds to 0.0625 : halfway between coil and head center
yview = round(Nvox/2);
xview = round(Nvox/2);
FOV = 0.5; % m

make_masks


current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

R1 = 0.06;  % m
aperture = pi/6;



% Calculate E field from coil alone:
[E1 Ex1 Ey1 Ez1] =  make_fig8(current, R1, xgap, 0.125, theta);
% bottom shield
[E2 Ex2 Ey2 Ez2] = make_doubleC(current, Rbottom, Zbottom, 0, pi/6, 0);
% top shield
[E3 Ex3 Ey3 Ez3] = make_doubleC(current, Rtop, Ztop, 0, pi/8, 0);

Ex = sfactorBottom*Ex2 + sfactorTop*Ex3 + Ex1;
Ey = sfactorBottom*Ey2 + sfactorTop*Ey3 + Ey1;
Ez = sfactorBottom*Ez2 + sfactorTop*Ez3 + Ez1;
% Ex = sfactorBottom*Ex2  + Ex1;
% Ey = sfactorBottom*Ey2  + Ey1;
% Ez = sfactorBottom*Ez2  + Ez1;

E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);

% calculate measures of cost
targetN = sum(targetmask(:));

Ehead = E .* headmask;
Ecyl = E .* cylmask;
sharp = sum(Ecyl(:)) / sum(Ehead(:))
prec =  sum(E(:).*targetmask(:)) / sum(Ehead(:));

% figure out penetration within cylinder of interest:
% i.e.:  where is abs(E) = 0.5*abs(E_surface) ?
Ecyl_max = max(Ecyl(:));
active = zeros(size(Ecyl));
active(Ecyl > Ecyl_max/2) = 1;
inds = find(active);
[x, y, z ]= ind2sub([Nvox, Nvox, Nvox], inds);
penetration = min(z);

hold off

E(E==max(E(:)))= 0;
ov([], (1+headmask+cylmask)*5 + E *256/(max(E(:)) - min(E(:))) ,...
    xview, yview, zview,0);
set(gcf, 'Name', sprintf('shield= %f pen= %f, sharp= %f', sfactorBottom, penetration, sharp))
drawnow
cost = norm([1/sharp])

subplot(224)
plot(squeeze(E(xview, yview,:)))
% + penetration

return


