%here are the commands to read the data ft and phase correct

[np,ns,nv,pss,lro,thk,data,seqfil,array, TI]=ReadFID5;

fprof_1 =fftshift(fft(data),1);
eps = 0.00001;
pc = conj(fprof_1(:,9)) ./ (abs(fprof_1(:,9) + eps));  %this is the phase correction part
for nir = 1:9
fprof_1(:,nir) = fprof_1(:,nir) .* pc;
end
plot(real(reshape(fprof_1,128*9,1)))
figure
plot(TI,real(sum(fprof_1(35:55,:))),'*')

% Now do the fitting:
decay = sum(real(fprof_1(35:55,:)));
Mo_guess = 1.2*decay(end);
T1_guess = 2;
TR = 15;

LB = [100, 0.5];
UB = [10e8, 3];

guess0 = [Mo_guess; T1_guess];
guess = lsqnonlin('ir_lsq', ...
    guess0, LB, UB, ...
    [],...
    TI, ...
    TR, ...
    decay);

M0 = guess(1)
T1 = guess(2)

m = ir_lsq(guess,TI,TR);
hold on
plot(TI,m,'r')
hold off