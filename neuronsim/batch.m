
hh
count=1;
for Amp=linspace(10,50,10)
for w=linspace(0,5,20)
    w
    Iext=pulse(0:dt:T,2,10,-Amp,w);
    t=0; hh;
    Vmax(count)=max(abs(V)); count = count+1;
    drawnow; pause(0.1)
end
end

Vmax = reshape(Vmax, 20,10);
figure; imagesc(Vmax)