function H=rsource(rs,js,H)

%%%%% this subroutine propagates fields all over the domain
%only values from 2 to nx are used.
ro=[0 0 0];
global nx ny nz dx dy dz;
for k=1:nz
for j=1:ny
for i=1:nx        
            ro0=[dx*(i-1/2);dy*(j-1/2);dz*(k-1)];
            ro1=[dx*(i-1);dy*(j-1/2);dz*(k-1/2)];
            ro2=[dx*(i-1);dy*(j-1/2);dz*(k-1/2)];
            node0=globalnode(i,j,k,0);
            node1=globalnode(i,j,k,1);
            node2=globalnode(i,j,k,2);
            for int=1:length(rs(:,1))
            R=-[rs(int,1)-ro0(1) rs(int,2)-ro0(2) rs(int,3)-ro0(3)];
            if abs(norm(R))~=0
            H(node0)=H(node0)+(js(int,1)*R(2)-js(int,2)*R(1))/(4*pi*norm(R)^3)*dx*dy;
            end
            R=-[rs(int,1)-ro1(1) rs(int,2)-ro1(2) rs(int,3)-ro1(3)];
            if abs(norm(R))~=0
            H(node1)=H(node1)+(js(int,2)*R(3)-js(int,3)*R(2))/(4*pi*norm(R)^3)*dy*dz;
            end           
            R=-[rs(int,1)-ro2(1) rs(int,2)-ro2(2) rs(int,3)-ro2(3)];
            if abs(norm(R))~=0
            H(node2)=H(node2)+(js(int,3)*R(1)-js(int,1)*R(3))/(4*pi*norm(R)^3)*dx*dz;
            end
            end
end
end
end