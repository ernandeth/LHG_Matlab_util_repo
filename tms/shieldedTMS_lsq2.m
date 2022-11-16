function [cost, E] = shieldedTMS_lsq2(parms, consts)
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
% W = relative weight of penetratin penalty (default=1)
%

global FOV Nvox;
global hst
global E1 Ex1 Ey1 Ez1


parms

doDivergence=0;

% sfactorBottom = parms(1);
% sfactorTop = parms(2);
% theta = parms(3);
% Rtop = parms(4);
% Rbottom = parms(5);
% Zbottom = parms(6);
% top_Z = parms(7);

sfactorBottom = parms(1);
Rbottom = parms(2);
Zbottom=parms(3);


% Rbottom = 0.12;
% Zbottom = 0.135;
% theta = 0;

R1 = 0.03;  % m
xgap = -(Rbottom-R1);

zview = round(Nvox/2)+8; 
yview = round(Nvox/2);
xview = round(Nvox/2);

target_z = 0.0375;  % This is halfway between center and surface
target_z = round(target_z * Nvox/FOV + Nvox/2);

make_masks

current =1e3/1e-4; % this is actually dI/dt in Amps/sec. realistic for TMS


% % Calculate E field from coil alone:
% [E1 Ex1 Ey1 Ez1] =  make_fig8(current, 0.03, 0, 0.085, 0);
% save defaultFig8.mat E1 Ex1 Ey1 Ez1
% % instead of calculating it, load it from a file
% load defaultFig8.mat

% bottom shield
% [E2 Ex2 Ey2 Ez2] = make_doubleC(current, Rbottom, Zbottom, 0, pi/6, 0);
% top shield
% [E3 Ex3 Ey3 Ez3] = make_doubleC(current, Rtop, top_Z, 0, pi/8, 0);
% [E1 Ex1 Ey1 Ez1] =  make_fig8(current, 0.03, 0, 0.085, 0);
% save defaultFig8.mat E1 Ex1 Ey1 Ez1

% figure 8 shield
[E2 Ex2 Ey2 Ez2] =  make_fig8(current, Rbottom, xgap, Zbottom, 0);

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
E1cyl = E1.* cylmask;

Eheadmax = max(Ehead(:));
Ecylmax = max(Ecyl(:));
E1cylmax = max(E1cyl(:));

sharp = sum(Ecyl(:)) / sum(Ehead(:))
%sharp = sum(Ecyl(Ecyl>=Ecylmax/2)) / sum(Ehead(Ehead>=Ecylmax/2))


% figure out penetration within cylinder of interest:
% i.e.:  where is abs(E) = 0.5*abs(E_surface) ?
% center =ceil([Nvox, Nvox, Nvox]/2);
penetration = E(xview,yview, target_z) / Ecylmax;

% W=1;  % weight of penetration on the cost function
if ~isfield(hst,'W')
    W=1;
else
    W = hst.W;
end

cost = norm([1/sharp  W/penetration]);

% keep track of the history of the algorithm:
hst.penetration = [hst.penetration ; penetration];
hst.sharp = [hst.sharp; sharp];
hst.parms = [hst.parms; parms];
hst.cost = [hst.cost; cost];
hst.Ecylmax = [hst.Ecylmax; Ecylmax];
	

if 0
    hold off
% 	E(E==max(E(:)))= 0;
% 	E1(E1==max(E1(:)))= 0;
	
    % normalization to the max on surface:
	E = E/Ecylmax;
	E1norm = E1/E1cylmax;
	
	sfigure(1);
	ov([], headmask + E*256/(Ecylmax) ,...
	    xview, yview, zview,0);
	set(gcf, 'Name', sprintf('sfac= %0.3f pen= %0.3f, sharp= %0.3f', sfactorBottom, penetration, sharp))
	drawnow



    p = get(gcf,'Position');
    p(1) = p(1)+500;

    sfigure(2); set(gcf,'Position',p);

    subplot (3,1,1)
    plot(squeeze(E(:,yview, zview))), title('x axis')
    hold on
    plot(squeeze(E1norm(:,yview, zview)),'r'), title('x axis')    
	legend('Total','Fig-8 only')
    hold off

    subplot (3,1,2)
    plot(squeeze(E(xview,:,zview))), title('y axis')
    hold on
    plot(squeeze(E1norm(xview,:,zview)),'r'), title('y axis')
    hold off

    subplot (3,1,3)
    plot(squeeze(E(xview, yview,:))), title('z axis')
    hold on
    plot(squeeze(E1norm(xview, yview,:)),'r'), title('z axis')
	axis([0 33 0 2])
    hold off

	sfigure(1);
end
return



