function H=countourxy(rs,js,i,j,k)
global nx ny dx dy dz;
H=0;
         rijk=[dx;dy;dz].*[(i-1);(j-1);(k-1)];
         
         ri1jk=[dx;dy;dz].*[(i);(j-1);(k-1)];
         
         rij1k=[dx;dy;dz].*[(i-1);(j);(k-1)];
         
         ri1j1k=[dx;dy;dz].*[(i);(j);(k-1)];

	     dx1=dh(rijk,ri1jk,rs,js,1);
	     dx2=dh(rij1k,ri1j1k,rs,js,1);
	     dy1=dh(rijk,rij1k,rs,js,2);
         dy2=dh(ri1jk,ri1j1k,rs,js,2);
if i~=1 && i~=nx && j~=1 && j~=ny
   H=H+dx1+dy2-dx2-dy1; 
elseif i==1 && j==1
        H=H-dx2+dy2;  
elseif i==1 && j~=1 && j~=ny
	     H=H+dx1-dx2+dy2;                 
elseif i==1 && j==ny
	     H=H+dx1+dy2;         
elseif i~=1 && i~=nx && j==1
	     H=H-dx2-dy1+dy2;        
elseif i~=1 && i~=nx && j==ny
	     H=H+dx1-dy1+dy2;          
elseif i==nx && j==1
	     H=H-dx2-dy1;
elseif i==nx && j~=1 && j~=ny
	     H=H+dx1-dx2-dy1;
elseif i==nx && j==ny
	     H=H+dx1-dy1;
end
