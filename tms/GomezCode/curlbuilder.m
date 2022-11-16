function curlmatrix=curlbuilder()


global nx ny nz dx dy dz;
col(4*curlequation(nx,ny,nz,2)-4*(nx*ny+nz*nx+nz*ny))=0;
row(4*curlequation(nx,ny,nz,2)-4*(nx*ny+nz*nx+nz*ny))=0;
value(4*curlequation(nx,ny,nz,2)-4*(nx*ny+nz*nx+nz*ny))=0;
sparseindex=1;
for i=1:nx
    for j=1:ny
        for k=1:nz
            ip=i+1;
            jp=j+1;
            kp=k+1;
            
   %%%%%%%%%%write the facexy equation
           if k~=1
            eqnfacecol=curlequation(i,j,k,0);
            if j~=1
            %Exijk*dx
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,0);
            value(sparseindex)=dx;
            sparseindex=sparseindex+1;  
            
            end
            if j~=ny 
            %Exij+1k*(-dx)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,jp,k,0);
            value(sparseindex)=-dx;
            sparseindex=sparseindex+1;
   
            end

            if i~=1 
            %Eyijk*(-dy)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,1);
            value(sparseindex)=-dy;
            sparseindex=sparseindex+1;
            end
            if i~=nx 
            %Eyi+1jk*(dy)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(ip,j,k,1);
            value(sparseindex)=dy;
            sparseindex=sparseindex+1;
 
            end
              
          end
 %%%%%facexy equation done%%%%%%%%%%%%%%
          if i~=1 
            %%%%%%%%%%write the faceyz equation
            eqnfacecol=curlequation(i,j,k,1);
            if k~=1 
            %Eyijk*dy
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,1);
            value(sparseindex)=dy;
            sparseindex=sparseindex+1;
             
            end
            
            if k~=nz  
            %Eyijk+1*(-dy)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,kp,1);
            value(sparseindex)=-dy;
            sparseindex=sparseindex+1;
   
            end
            if j~=1 
            %Ezijk*(-dz)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,2);
            value(sparseindex)=-dz;
            sparseindex=sparseindex+1;
    
            end
            if j~=ny 
            %Ezij+1k*(dz)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,jp,k,2);
            value(sparseindex)=dz;
            sparseindex=sparseindex+1;
 
            end
           end
      %%%%%%%%%%write the facexz equation
        if j~=1
            eqnfacecol=curlequation(i,j,k,2);
            if i~=1 
            %Ezijk*dz
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,2);
            value(sparseindex)=dz;
            sparseindex=sparseindex+1;
    
            end
            
            if i~=nx 
            %Ezi+1jk*(-dz)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(ip,j,k,2);
            value(sparseindex)=-dz;
            sparseindex=sparseindex+1;
   
            end
            if k~=1 
            %Exijk*(-dx)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,k,0);
            value(sparseindex)=-dx;
            sparseindex=sparseindex+1;
            
            end
            
            if k~=nz
            %Exijk+1*(dx)
            col(sparseindex)=eqnfacecol;
            row(sparseindex)=efieldnode(i,j,kp,0);
            value(sparseindex)=dx;
            sparseindex=sparseindex+1;
   
           end
      end
        end
    end
end
curlmatrix=sparse(col,row,value);



