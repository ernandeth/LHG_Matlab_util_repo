function n=efieldnode(i,j,k,dim)

global nx ny nz;
if dim==0
    n=i+(j-2)*nx+(k-2)*nx*(ny-1);
elseif dim==1
    n=i-1+(j-1)*(nx-1)+(k-2)*(nx-1)*ny +nx*ny*nz-nx*ny-nx*nz+nx;
elseif dim==2
    n=i-1+(j-2)*(nx-1)+(k-1)*(nx-1)*(ny-1)+2*nx*ny*nz-nx*ny-nx*nz-ny*nz-ny*nx+nx+ny;
end