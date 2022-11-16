beta = rand(4);
x = randn(20,4);

y=zeros(size(x));

for n=1:4
    y(:,n) = x * beta(:,n);
end

plot(y)
