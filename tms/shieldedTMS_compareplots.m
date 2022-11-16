close all
clear all


FOV=0.24;
Nvox=128;
zview = round(Nvox/2 + 0.0325*Nvox/FOV); % corresponds to 0.0625 : halfway between coil and head center
yview = round(Nvox/2);
xview = round(Nvox/2);

zmax= round(Nvox/2 + 0.065*Nvox/FOV);

load defaultFig8_128.mat
E0 = 20*log10(E1 /E1(xview, yview, zmax));

load (['fig8_shield_128_W10.mat'])
E1 = 20*log10(E2/E2(xview, yview, zmax));


load (['fig8_shield_128_W0.1.mat'])
E2 = 20*log10(E2/E2(xview, yview, zmax));




cmax = 0;
cmin = 20*log10(0.5);
cmin = 20*log10(1e-2);

sfigure(3);
ov02([],  E0 ,...
    xview, yview, zview,0, cmin, cmax);
set(gcf, 'Name',['No shield'])



sfigure(4);
ov02([],  (E1) ,...
    xview, yview, zview,0, cmin, cmax);


sfigure(5);
ov02([],  (E2) ,...
    xview, yview, zview,0, cmin, cmax);

