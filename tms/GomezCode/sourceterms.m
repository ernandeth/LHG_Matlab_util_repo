function Hsource=sourceterms(rs,js,H)
%Hsource=sourceterms(H) returns a vector Hsource, with the same length
%as H, but where the first nx entries
%have Hz the next ny entries have nx and the last entries have Hy
%used as source of the curl equations
mu0=1.25663706*10^-6;
mu0=1;
global nx ny nz dx dy dz;
tic
Hsource(1:curlequation(nx,ny,nz,2))=0;
for i=1:nx
    for j=1:ny
        for k=1:nz
if k~=1
    if i==1 || j==1 || i==nx || j==ny
    Hsource(curlequation(i,j,k,0))=(-mu0/(4*pi))*countourxy(rs,js,i,j,k);
    else
          Hsource(curlequation(i,j,k,0))=-mu0*H(globalnode(i,j,k,0));
    end
end
if i~=1
        if k==1 || j==1 || k==nz || j==ny
    Hsource(curlequation(i,j,k,1))=(-mu0/(4*pi))*countouryz(rs,js,i,j,k);
    else

        Hsource(curlequation(i,j,k,1))=-mu0*H(globalnode(i,j,k,1));
        end
end
if j~=1
    if i==1 || k==1 || i==nx || k==nz
    Hsource(curlequation(i,j,k,2))=(-mu0/(4*pi))*countourxz(rs,js,i,j,k);
    else
        Hsource(curlequation(i,j,k,2))=-mu0*H(globalnode(i,j,k,2));
    end
end



      end
    end
    end
toc