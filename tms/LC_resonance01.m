c = 10.^[-8:0.1:-2];
L1 = 150e-6;
f = 1./(2*pi*sqrt(L1*c));
plot(c,f)
xlabel('Cap. (F)')
ylabel('Res. Freq (Hz)')
axis([0 1e-3 0 5e3])



% checking out cylindrical inductors:
r = [0.01:0.001:0.05];
N = 20;
d = 1e-2;
L = (r.^2 * N^2) *1e-5 ./ (2.*r + 2.8*d)
plot(r,L)
