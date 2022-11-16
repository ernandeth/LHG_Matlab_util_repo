
n = zeros(1000,10);
for iter=1:1000
    
    nth=1;
    for th=linspace(0,1,10)
        a = randn(1000,1);
        a(a<th) = 0;
        a(a>th) = 1;
        
        b = randn(1000,1);
        b(b<th) = 0;
        b(b>th) = 1;

        b = b & a;
        
        nth = nth+1;
        
        n(iter, nth) = ncsty(a,b);
        
    end
end
n = sum(n,2)/nth;
hist(n,50)
