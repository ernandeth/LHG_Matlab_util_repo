% this function will simulate the spiral figure eight coil design
[x,y] = make_spiral_coil(2,1.1,6);
 x1 = x'/10+0.55; y1 = y'/10;

[x,y] = make_spiral_coil(2,1.1,6);
x2 = -x1; y2 = y1;
 
subplot(221)
plot(x1,y1); hold on; plot(x2,y2,'r'); hold off
axis equal
drawnow

wire1 = [x1,y1, zeros(size(x1))];

[E1, Ex1, Ey1 , Ez1]=Efield(1, wire1, [64 64 64], 1.5);
subplot(224)
lightbox(E1);
drawnow

wire2 = [x2,y2, zeros(size(x2))];

[E2, Ex2, Ey2 , Ez2]=Efield(1, wire2, [64 64 64], 1.5);
subplot(223)
lightbox(E2);
drawnow

Ext = Ex1 + Ex2;
Eyt = Ey1 + Ey2;
Ezt = Ez1 + Ez2;

Et = sqrt(Ext.^2 + Eyt.^2 + Ezt.^2);
subplot(222)
lightbox(Et);
drawnow


return
[x,y] = make_spiral_coil(2,1.1,12);

plot(x,y,'k'); fatlines(3); axis off ; axis square
print -dpng spiral

subplot(224)
plot(x,y,'k'); ; axis off ; axis square


subplot(221)
plot(-x,y,'k'); ; axis off ; axis square
subplot(223)
plot(-x,y,'k'); ; axis off ; axis square


