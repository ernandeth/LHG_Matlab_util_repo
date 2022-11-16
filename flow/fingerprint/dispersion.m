
for P = [0:0.3:4]
Disp = 10.00
%P = 1.5
transit = 0;
dt = 0.01;

ttmp = linspace(0,5, 5/dt)';
art_kernel = (ttmp-transit).^(Disp*P) .*  exp(-Disp *(ttmp-transit));
art_kernel(art_kernel<0)=0;
art_kernel = art_kernel / sum(art_kernel);  %normalize the mass in the dispersion

plot(art_kernel)
hold on
end