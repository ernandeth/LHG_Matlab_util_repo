clear all;
pecfolder='/home/luisgo/Finitedifferenceintegralwall/';
boxndim=[20;20;20];
boxddim=[.0017;.0017;.003];

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
r=.032;
R=.25;
zsource=-1;
apperture=pi/6;
distance=20;
[rs{1} Jvec{1}]=loopsource([boxndim(1)*boxddim(1)/2 boxndim(2)*boxddim(2)/2+r/2 (zsource)*boxddim(3)],[0 0 1],r,32);
[rs{2} Jvec{2}]=loopsource([boxndim(1)*boxddim(1)/2 boxndim(2)*boxddim(2)/2-r/2 (zsource)*boxddim(3)],[0 0 1],r,32);
Jvec{2}=Jvec{2}*(-1);
    wx = R*sin(apperture:2*pi/32:(2*pi-apperture));
    wy = -R*cos(apperture:2*pi/32:(2*pi-apperture));
    wz = dz*(zsource+distance)+wy*sin(0);   

rs{3}=[[wx(1) R+wy(1) -100*dz];
   [wx' R+wy' wz'];
    [wx(end) R+wy(end) -100*dz]];
rs{4}=[[wx(1) -R-wy(1) -100*dz];
    [wx' -R-wy' wz'];
    [wx(end) -R-wy(end) -100*dz]];
for j=3:4
    slen=length(rs{j});
    totlen=0;
rs{j}(:,1)=rs{j}(:,1)+nx/2*dx;
rs{j}(:,2)=rs{j}(:,2)+ny/2*dy;
for i=1:slen
Jvec{j}(i,:)=(rs{j}(i*(i~=slen)+1,:)-rs{j}(i,:));
    totlen=totlen+norm(Jvec{j}(i,:));
end
Jvec{j}=.25*Jvec{j}/totlen;
end
Jvec{4}=Jvec{4}*(-1);
%seesetup(rs,sigmap,ones(size(sigmap)))
Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec);


save(strcat(pecfolder,'runfinitecode.mat'))
