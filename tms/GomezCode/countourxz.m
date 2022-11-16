function H=countourxz(rs,js,i,j,k)
global nx ny nz dx dy dz;
H=0;

            rijk=[dx;dy;dz].*[(i-1);(j-1);(k-1)];
            ri1jk=[dx;dy;dz].*[(i);(j-1);(k-1)];
            rijk1=[dx;dy;dz].*[(i-1);(j-1);(k)];
            ri1jk1=[dx;dy;dz].*[(i);(j-1);(k)];
            
  
         dx1=dh(rijk,ri1jk,rs,js,1);
	     dx2=dh(rijk1,ri1jk1,rs,js,1);
	     dz1=dh(rijk,rijk1,rs,js,3);
         dz2=dh(ri1jk,ri1jk1,rs,js,3);
if k~=1 && k~=nz && i~=1 && i~=nx
   H=H+dz1+dx2-dz2-dx1; 
elseif k==1 && i==1
        H=H-dz2+dx2;  
elseif k==1 && i~=1 && i~=nx
	     H=H+dz1-dz2+dx2;                 
elseif k==1 && i==nx
	     H=H+dz1+dx2;         
elseif k~=1 && k~=nz && i==1
	     H=H-dz2-dx1+dx2;        
elseif k~=1 && k~=nz && i==nx
	     H=H+dz1-dx1+dx2;          
elseif k==nz && i==1
	     H=H-dz2-dx1;
elseif k==nz && i~=1 && i~=nx
	     H=H+dz1-dz2-dx1;
elseif k==nz && i==nx
	     H=H+dz1-dx1;
end
