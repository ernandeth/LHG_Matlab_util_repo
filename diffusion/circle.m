function result = circle(r)

w = 0:pi/100:2*pi;
x = r*cos(w);
y = r*sin(w);

result = plot(x,y);
return