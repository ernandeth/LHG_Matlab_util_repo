function z = tilen(data,dim);
% function z = tilen(data,dim);
% data should be 3- or 4-D, with dimensions a*b*c*d
% output is a 2d array with a*b images tiled in a dim(1)*dim(2) grid

n1 = size(data,1);
n2 = size(data,2);
n3 = size(data,3);
n4 = size(data,4);

if ~exist('dim','var')
    dim(1) = n3;%max(n3,n4);
    dim(2) = n4;%min(n3,n4);
end
n = n3*n4;
data = reshape(data,[n1,n2,n]);
z=[];
zz=[];
for k=1:n
    if mod(k,dim(1))==1 || dim(1)==1
        zz=[];
    end
    zz = cat(2,zz,data(:,:,k));
    
    if k==n && mod(k,dim(1))
        zz = cat(2,zz,zeros(n1,n2*(dim(1)-mod(k,dim(1)))));
    end
        
    if mod(k,dim(1))==0 || k==n
        z = cat(1,z,zz);
    end
end