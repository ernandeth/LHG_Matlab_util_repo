function Efield=oneshotcerry(boxndim,boxddim,sigmap,rs,Jvec)

global nx ny nz dx dy dz;
nx=boxndim(1);
ny=boxndim(2);
nz=boxndim(3);
dx=boxddim(1);
dy=boxddim(2);
dz=boxddim(3);
% plot3(rs{1}(:,1),rs{1}(:,2),rs{1}(:,3))
% figure
H=zeros([globalnode(nx,ny,nz,2) 1]);

% generate H field in the domain


%%% generate source vector including special boundary terms
tic
Hsource=contoursourceterms(rs,Jvec);
toc
sigmaav=sigmaaverage(sigmap);
cond=condmatrix(sigmaav);
curl=curlbuilder();
sigmasol=zeros(length(cond(:,1)),1);



%create system matrix and source terms
A=cat(1,cond,curl);
b=cat(1,sigmasol, transpose(Hsource));
tic
x=lsqr(A,b,10^-7,2000); %solve
toc
%extract solution
en=0;
st=1;
en=nx*(ny-1)*(nz-1);
Efield{1}=reshape(x(st:en),[nx ny-1 nz-1]);
st=1+en;
en=en+ny*(nx-1)*(nz-1);
Efield{2}=reshape(x(st:en),[nx-1 ny nz-1]);
st=1+en;
en=en+nz*(nx-1)*(ny-1);
Efield{3}=reshape(x(st:en),[nx-1 ny-1 nz]);
