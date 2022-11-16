function cond=condmatrix(sigmaavg)
global nx ny nz dx dy dz;

col(6*(nx-1)*(ny-1)*(nz-1))=0;
row(6*(nx-1)*(ny-1)*(nz-1))=0;
value(6*(nx-1)*(ny-1)*(nz-1))=0;

sparseindex=1;
eqnfacecol=1;
for k=2:nz
    for j=2:ny
        for i=2:nx
            ip=i-1;
            jp=j-1;
            kp=k-1;
            

            %%%%%%%%%%write the Ex part of equation

            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(ip,j,k,0);
            value(sparseindex)=-sigmaavg(ip,j,k,1)*dy*dz;
            sparseindex=sparseindex+1;
            %-sigmaijk*Exi+1jk
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,0);
            value(sparseindex)=sigmaavg(i,j,k,1)*dy*dz;
            sparseindex=sparseindex+1;
 
            %%%%%%%%%%write the Ey part of equation

            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,jp,k,1);
            value(sparseindex)=-sigmaavg(i,jp,k,2)*dx*dz;
            sparseindex=sparseindex+1;
            %-sigmaij+1k*Eyijk
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,1);
            value(sparseindex)=sigmaavg(i,j,k,2)*dx*dz;
            sparseindex=sparseindex+1;
            %%%%%%%%%%write the Ez part of equation

            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,kp,2);
            value(sparseindex)=-sigmaavg(i,j,kp,3)*dx*dy;
            sparseindex=sparseindex+1;
            %-sigmaijk*Ezijk
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,2);
            value(sparseindex)=sigmaavg(i,j,k,3)*dx*dy;
            sparseindex=sparseindex+1;
           

             eqnfacecol=eqnfacecol+1;

       
        end
    end
end

cond=sparse(col,row,value,eqnfacecol-1,efieldnode(nx,ny,nz,2)); %ensure dimensions are right create matrix