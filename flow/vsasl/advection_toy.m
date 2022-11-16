Ma = zeros(100,1);
Ma(2) = 1;
for n = 2:1000;
    
    for m =2:100
        dMa = 10*(Ma(m-1)-Ma(m));
        Ma(m) = Ma(m) + dMa*dt;
   
    end
         plot(Ma,'r')
         axis([1 10 0 1])
        drawnow
end