function sim_active_shield(R1, R2, L1, L2, N1, N2, F1, F2)
%%
close all

R1 = 0.5; R2 = 1.5*R1;
L1 = 0.5; L2 = 0.15*L1;

N1 = 10; N2 = 5;
F1 = 1; F2 = -0.3;

[x y z] = make_solenoid(R1, L1, N1, 200);
% Two coils at the ends of the main coil
[x2 y2 z2] = make_solenoid(R2, L2, N2, 200);
z2 = z2+L1/2;
z3 = z2-2*L1/2;
hold off



[b1 b1x b1y b1z] = biot3d(F1, [x y z]);

[b2 b2x b2y b2z] = biot3d(F2, [x2 y2 z2]);

b1(isnan(b1))=0;
b2(isnan(b2))=0;

b3 = b2(:, :, end:-1:1);
b3x = b2x(:, :, end:-1:1);
b3y = b2y(:, :, end:-1:1);
b3z = b2z(:, :, end:-1:1);

%%
% Normalize them:
Bscale = abs(b1(round(end/2), round(end/2), round(end/2)));
b1x = b1x/Bscale;
b1y = b1y/Bscale;
b1z = b1z/Bscale;
b1 = b1/Bscale;

Bscale = abs(b2(round(end/2), round(end/2), round(end/2)));
b2x = b2x/Bscale;
b2y = b2y/Bscale;
b2z = b2z/Bscale;
b2 = b2/Bscale;

Bscale = abs(b3(round(end/2), round(end/2), round(end/2)));
b3x = b3x/Bscale;
b3y = b3y/Bscale;
b3z = b3z/Bscale;
b3 = b3/Bscale;

SF = 0.15;
close all
figure
subplot(211)
% imagesc((abs(squeeze(0.5*b1z(100,:, :) ))));
% caxis([-1 1]*2) 
% title('Reduced Main')
%caxis auto
plot3(x, z, y, 'k')
hold on
plot3(x2, z2, y2,'r')
plot3(x2, z3, y2, 'r')
axis square
axis off

subplot(234)
imagesc(((squeeze(-b1z(100,:, :)))));
caxis([-1 1]*2) 
title('Main Field')
axis square
%caxis auto
axis off

subplot(235)
imagesc(((squeeze(SF*(b3z(100,:,:) + b2z(100,:, :))))));
caxis([-1 1]*2) 
axis square
title('Shielding Field')
%caxis auto
axis off


subplot(236)
imagesc(((squeeze(-b1z(100,:, :) + SF*(b2z(100,:, :) + b3z(100,:, :))))));
caxis([-1 1]*2) 
axis square
title('Net Field')
%caxis auto
axis off

colormap jet