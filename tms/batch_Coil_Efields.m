
zview = Nvox/2 +8; % corresponds to 0.0625 : halfway between coil and head center
yview = Nvox/2;
xview = Nvox/2;


current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

config='hemisphere'
r = 0.06;  % this radius is a fraction of the FOV
aperture = pi/6;

R = 0.10;
shield_offset = .0625;
theta = pi/8;

FOV = 0.5; % m
% this means that each voxel's size  is 0.5/64 meters = 0.0078 m

 %[E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, theta);
    [E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, 0);
    save Efield_fig8_01 E1 Ex1 Ey1 Ez1
    [E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, pi/4);
    save Efield_fig8_02 E1 Ex1 Ey1 Ez1
    [E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, pi/6);
    save Efield_fig8_03 E1 Ex1 Ey1 Ez1
    [E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, pi/8);
    save Efield_fig8_04 E1 Ex1 Ey1 Ez1
    [E1 Ex1 Ey1 Ez1] =  make_fig8(current, r, 0.125, -pi/8);
    save Efield_fig8_05 E1 Ex1 Ey1 Ez1
    
    % big lower shield is #2
    %[E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.07, 0.125 - shield_offset, 0, pi/6, 0);
    [E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.07, 0.125 - shield_offset/2, 0, pi/6, 0);
    save Efield_bottomshield_01 E2 Ex2 Ey2 Ez2
    [E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.07, 0.125 - shield_offset, 0, pi/6, 0);
    save Efield_bottomshield_02 E2 Ex2 Ey2 Ez2
    [E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.10, 0.125 - shield_offset/2, 0, pi/6, 0);
    save Efield_bottomshield_03 E2 Ex2 Ey2 Ez2
    [E2 Ex2 Ey2 Ez2] = make_doubleC(current, 0.06, 0.125 - shield_offset/2, 0, pi/6, 0);
    save Efield_bottomshield_04 E2 Ex2 Ey2 Ez2
    
    % little top shield is #3
    %[E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/8, 0);
    [E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/8, 0);
    save Efield_topshield_01 E3 Ex3 Ey3 Ez3
    [E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.03, 0.125 + shield_offset/2, 0, pi/6, 0);
    save Efield_topshield_02 E3 Ex3 Ey3 Ez3
    [E3 Ex3 Ey3 Ez3] = make_doubleC(current, 0.05, 0.125 + shield_offset/2, 0, pi/8, 0);
    save Efield_topshield_03 E3 Ex3 Ey3 Ez3
    
