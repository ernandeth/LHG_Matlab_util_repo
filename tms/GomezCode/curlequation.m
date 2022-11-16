function n=curlequation(i,j,k,dim)
%returns index value for the i,j,k,dim number
global nx ny nz;
if dim==0
n=i+(j-1)*nx+(k-2)*nx*ny;
elseif dim==1
    n=i-1+(j-1)*(nx-1)+(k-1)*(nx-1)*ny+dim*nx*ny*nz-nx*ny;
elseif dim==2
    n=i+(j-2)*nx+(k-1)*nx*(ny-1)+dim*nx*ny*nz-nx*ny-ny*nz;
end