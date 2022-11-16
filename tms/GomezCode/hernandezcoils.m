clear all;
pecfolder='/home/luisgo/Finitedifferenceintegralwall/';
boxndim=[20;20;20];
boxddim=[.0019;.0019;.003];

sigmap=sigma3(1);
brainlx=size(sigmap);
p=6;
sigma(1:p+brainlx(1),1:p+brainlx(2),1:p+brainlx(3))=0;
sigma(p/2+1:p/2+brainlx(1),p/2+1:p/2+brainlx(2),p/2+1:p/2+brainlx(3))=sigmap;
sigmap=sigma;
clear sigma;
boxndim=size(sigmap);
global nx ny nz dx dy dz;
nx=boxndim(1);
ny=boxndim(2);
nz=boxndim(3);
dx=boxddim(1);
dy=boxddim(2);
dz=boxddim(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0.06;  % this radius is a fraction of the FOV
aperture = pi/6;
zpos=p/2*dz-.010;
R = 0.10;
sfactor = -0.33;
shield_offset = .0625;
theta=pi/8;

current =1e3/1e-4; % this is actually dI/dt in Amps/sec.

%%%%%%%%%%%%%%%%%first coil simulation
[rs{1},Jvec{1}] = make_fig8(current, r, zpos, theta)
for i=1:length(rs)
    rs{i}(:,1)=rs{i}(:,1)+nx*dx/2;
    rs{i}(:,2)=rs{i}(:,2)+ny*dy/2;
end
Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec);
[Ex1 Ey1 Ez1]=organizefield(Efield);
save(strcat(pecfolder,'Efield1.mat'),'Ex1','Ey1','Ez1')

%%%%%%%%%%%%%%%%%second coil simulation
[rs{1},Jvec{1}] = make_doubleC(current, 0.087, zpos + shield_offset, 0, pi/6);
for i=1:length(rs)
    rs{i}(:,1)=rs{i}(:,1)+nx*dx/2;
    rs{i}(:,2)=rs{i}(:,2)+ny*dy/2;
end

Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec);
[Ex2 Ey2 Ez2]=organizefield(Efield);
save(strcat(pecfolder,'Efield2.mat'),'Ex2','Ey2','Ez2')

%%%%%%%%%%%%%%%%%third coil simulation
[rs{1},Jvec{1}] = make_doubleC(current, 0.03, zpos - shield_offset/2, 0, pi/8);
for i=1:length(rs)
    rs{i}(:,1)=rs{i}(:,1)+nx*dx/2;
    rs{i}(:,2)=rs{i}(:,2)+ny*dy/2;
end
Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec);
[Ex3 Ey3 Ez3]=organizefield(Efield);
save(strcat(pecfolder,'Efield3.mat'))
