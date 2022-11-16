function n=globalnode(i,j,k,dim)
%returns index value for the i,j,k,dim number
global nx ny nz;

n=i+(j-1)*nx+(k-1)*nx*ny+dim*nx*ny*nz;