function H=countouryz(rs,js,i,j,k)
global nx ny nz dx dy dz;
H=0;
            rijk=[dx;dy;dz].*[(i-1);(j-1);(k-1)];
            rij1k=[dx;dy;dz].*[(i-1);(j);(k-1)];
            rijk1=[dx;dy;dz].*[(i-1);(j-1);(k)];
            rij1k1=[dx;dy;dz].*[(i-1);(j);(k)];
         dy1=dh(rijk,rij1k,rs,js,2);
	     dy2=dh(rijk1,rij1k1,rs,js,2);
	     dz1=dh(rijk,rijk1,rs,js,3);
         dz2=dh(rij1k,rij1k1,rs,js,3);
    
if k~=1 && k~=nz && j~=1 && j~=ny
   H=H+dy1+dz2-dy2-dz1; 
elseif j==1 && k==1
        H=H-dy2+dz2;  
elseif j==1 && k~=1 && k~=nz
	     H=H+dy1-dy2+dz2;                 
elseif j==1 && k==nz
	     H=H+dy1+dz2;         
elseif j~=1 && j~=ny && k==1
	     H=H-dy2-dz1+dz2;        
elseif j~=1 && j~=ny && k==nz
	     H=H+dy1-dz1+dz2;          
elseif j==ny && k==1
	     H=H-dy2-dz1;
elseif j==ny && k~=1 && k~=nz
	     H=H+dy1-dy2-dz1;
elseif j==ny && k==nz
	     H=H+dy1-dz1;
end
