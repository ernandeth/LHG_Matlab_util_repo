function [Ex Ey Ez]=organizefield(Efield);
global nx ny nz;

Ex(1:nx,1,1:nz)=zeros([nx 1 nz]);
Ex(1:nx,1:ny,1)=zeros([nx ny 1]);
Ex(1:nx,2:ny,2:nz)=Efield{1};

Ey(1,1:ny,1:nz)=zeros([1 ny nz]);
Ey(1:nx,1:ny,1)=zeros([nx ny 1]);
Ey(2:nx,1:ny,2:nz)=Efield{2};

Ez(1:nx,1,1:nz)=zeros([nx 1 nz]);
Ez(1,1:ny,1:nz)=zeros([1 ny nz]);
Ez(2:nx,2:ny,1:nz)=Efield{3};

