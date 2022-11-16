function [cost, E] = shieldedTMS_lsq3(parms, consts)
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
% this version used a square shield
global FOV Nvox
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
Zbottom = parms(3);

% FOV = 0.5; % m
% Nvox=33;
zview = round(Nvox/2 + Nvox/8); % corresponds to 0.0625 : halfway between coil and head center
yview = round(Nvox/2);
xview = round(Nvox/2);

make_masks


current =1e3/1e-4; % this is actually dI/dt in Amps/sec.
% some hard coded parms:
xgap = 0;
R1 = 0.06;  % m
aperture = pi/6;
theta = 0;
%Zbottom = 0.135;


% Calculate E field from coil alone:
 [E1 Ex1 Ey1 Ez1] =  make_fig8(current, R1, xgap, 0.125, theta);
%load defaultFig8.mat

% bottom shield
[E2 Ex2 Ey2 Ez2] = make_square_shield(current, Rbottom, Zbottom);

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
center =ceil([Nvox, Nvox, Nvox]/2);
penetration = E(center(1),center(2), center(3)) / Ecylmax

cost = norm([1/sharp  1/penetration])
%cost = 1/sharp;
	
hold off
if 1

% 	E(E==max(E(:)))= 0;
% 	E1(E1==max(E1(:)))= 0;
	
	E = E/Ecylmax;
	E1 = E1/E1cylmax;
	
	sfigure(1);
	ov([], (1+headmask+cylmask)*5 + E *128/(max(E(:)) - min(E(:))) ,...
	    xview, yview, zview,0);
	set(gcf, 'Name', sprintf('shield= %f pen= %f, sharp= %f', sfactorBottom, penetration, sharp))
	drawnow



    p = get(gcf,'Position');
    p(1) = p(1)+500;

    sfigure(2); set(gcf,'Position',p);

    subplot (3,1,1)
    plot(squeeze(E(:,yview, zview))), title('x axis')
    hold on
    plot(squeeze(E1(:,yview, zview)),'r'), title('x axis')    
	legend('Total','Fig-8 only')
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
	axis([0 33 0 2])
    hold off

	sfigure(1);
end
return


