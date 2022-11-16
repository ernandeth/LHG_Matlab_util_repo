% this function will simulate the spiral figure eight coil design
[x,y] = make_spiral_coil(2,1.1, 15);


plot(x,y,'k'); 


axis square
axis([-6.8 9.8 -8.1 8.4])
axis off

fatlines(3)
drawnow

print -deps spiral_15turn


