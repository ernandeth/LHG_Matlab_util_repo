N=10;
y=randn(N,1);
F=randn(N,N);
B=randn(N,1);
s=randn(N,1);
M=randn(N,1);

J1=0.5*norm(y-(F.*exp(-j*B*s'))*M)^2

J2=0.5*y'*y-0.5*y'*((F.*exp(-j*B*s'))*M)-0.5*((F.*exp(-j*B*s'))*M)'*y+...
    0.5*((F.*exp(-j*B*s'))*M)'*((F.*exp(-j*B*s'))*M)

tmp=0;
tmp1=0;
tmp2=0;

for i=1:N,
    for k=1:N
        tmp = tmp + 0.5*y(i)^2/N;
        tmp1 = tmp1 + y(i)*F(i,k)*cos(B(i)*s(k))*M(k);
    end
end

for i=1:N,
    for k=1:N
        for kd=1:N
            tmp2 = tmp2 + 0.5*F(i,k)*F(i,kd)*( cos(B(i)*(s(k)-s(kd))) )*M(k)*M(kd);
        end
    end
end

J3 = tmp - tmp1 + tmp2
