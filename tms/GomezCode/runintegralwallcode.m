clear all
global nx ny nz dx dy dz;
nx=10;
ny=10;
nz=10;
dx=.001;
dy=.001;
dz=.001;
sigmap=zeros([nx ny nz]);
[rs{1} Jvec{1}]=loopsource([nx*dx/2 ny*dy/2 2*nz*dz],[0 0 1],.010,32);
% plot3(rs{1}(:,1),rs{1}(:,2),rs{1}(:,3))
% figure
H=zeros([globalnode(nx,ny,nz,2) 1]);

for i=1:length(rs)
H=rsource(rs{i},Jvec{i},H);
end

% H=reshape(H,[nx ny nz 3]);
% subplot(1,3,1),imagesc(H(:,:,5,1));
% subplot(1,3,2),imagesc(H(:,:,5,2));
% subplot(1,3,3),imagesc(H(:,:,5,3));
% H=H(:);

hvect=sourceterms(rs,Jvec,H);

figure
en=0;
st=1;
en=nx*ny*(nz-1);
hsource{3}=reshape(hvect(st:en),[nx ny nz-1]);
st=1+en;
en=en+nz*ny*(nx-1);
hsource{1}=reshape(hvect(st:en),[nx-1 ny nz]);
st=1+en;
en=en+nz*nx*(ny-1);
hsource{2}=reshape(hvect(st:en),[nx ny-1 nz]);

subplot(1,3,1),imagesc(hsource{1}(:,:,5));
subplot(1,3,2),imagesc(hsource{2}(:,:,5));
subplot(1,3,3),imagesc(hsource{3}(:,:,5));
% 
% 
% sigmaav=sigmaaverage(sigmap);
% cond=condmatrix(sigmaav);
% curl=curlbuilder();
% sigmasol=zeros(length(cond(:,1)),1);
% 
% 
% 
% create system matrix and source terms
% A=cat(1,cond,curl);
% b=cat(1,sigmasol, transpose(hvect));
% 
% x=lsqr(A,b,10^-7,500);
% 
% 
% en=0;
% st=1;
% en=nx*(ny-1)*(nz-1);
% Efield{1}=reshape(x(st:en),[nx ny-1 nz-1]);
% st=1+en;
% en=en+ny*(nx-1)*(nz-1);
% Efield{2}=reshape(x(st:en),[nx-1 ny nz-1]);
% st=1+en;
% en=en+nz*(nx-1)*(ny-1);
% Efield{3}=reshape(x(st:en),[nx-1 ny-1 nz]);
% 
% imagesc(Efield{3}(:,:,5))