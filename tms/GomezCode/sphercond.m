function sigma=sphercond(rcenter,a1,a2,a3,cond1,cond2,cond3,sigma)
%returns sigma profile in mesh rlat of a cylinder object with conductivity=cond centered
%about robj and radius=a
global nx ny nz dx dy dz;


noderadiusp=@(m,n,p) [(m-1/2)*dx (n-1/2)*dy (p-1/2)*dz];


for i=1:nx
    for j=1:ny
        for k=1:nz
        nodelocation=noderadiusp(i,j,k);
                if sqrt((nodelocation(1)-rcenter(1))^2+(nodelocation(2)-rcenter(2))^2+(nodelocation(3)-rcenter(3))^2)<=a1   %if inside cylinder assign loss of 1S/m
                sigma(i,j,k)=sigma(i,j,k)+cond1;
                elseif sqrt((nodelocation(1)-rcenter(1))^2+(nodelocation(2)-rcenter(2))^2+(nodelocation(3)-rcenter(3))^2)<=a2
                sigma(i,j,k)=sigma(i,j,k)+cond2;    
                elseif sqrt((nodelocation(1)-rcenter(1))^2+(nodelocation(2)-rcenter(2))^2+(nodelocation(3)-rcenter(3))^2)<=a3
                sigma(i,j,k)=sigma(i,j,k)+cond3;       
               end
        end
    end
end

