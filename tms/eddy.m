t = [1:1000]/1000;
R = 40;
a = zeros(1,1000);
a(101:150) = (1:50)/50;
a(151:800) = 1;        
a(801:850) = (50:-1:1)/50;
g = a;
plot(g)
hold on
e = exp(-t*R);
ge = conv(g,e);
plot(ge)
c = exp(t*R)./(2*t);
plot(c)
gec = conv(ge,c);
plot(gec,'k')