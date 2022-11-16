function sigmaavg=sigmaaverage(sigma)
%takes average of the conductivities at each node
%sigmaav=1/4[sigma(i-1,j,k)+sigma(i,j-1,k)+sigma(i,j,k-1)+sigma(i,j,k)]
global nx ny nz;
sigmaavg(nx,ny,nz,3)=0;
for i=1:nx
    for j=1:ny
        for k=1:nz
            if i==1
                ip=i;
            else
                ip=i-1;
            end
            if j==1
                jp=j;
            else
                jp=j-1;
            end
            if k==1
                kp=k;
            else
                kp=k-1;
            end
     
           
sigmaavg(i,j,k,1)=(1/4)*(sigma(i,j,k)+sigma(i,jp,k)+sigma(i,j,kp)+sigma(i,jp,kp));
sigmaavg(i,j,k,2)=(1/4)*(sigma(i,j,k)+sigma(ip,j,k)+sigma(i,j,kp)+sigma(ip,j,kp));
sigmaavg(i,j,k,3)=(1/4)*(sigma(i,j,k)+sigma(ip,j,k)+sigma(i,jp,k)+sigma(ip,jp,k));

        end
    end
end
