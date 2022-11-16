%read the t2 profile stuff
%bm july 2005
[np,ns,nv,pss,lro,thk,data,seqfil,array]=ReadFID5;
TE = [15:15:32*15]/1000;
fprof_1 =fftshift(fft(data),1);
eps = 0.00001;

%phase correct
pc = conj(fprof_1(:,2)) ./ (abs(fprof_1(:,2) + eps));
for nir = 1:32
fprof_1(:,nir) = fprof_1(:,nir) .* pc;
end

plot(abs(reshape(fprof_1,128*32,1)));

array_val = 15*[1:32];

figure

plot(TE(2:2:end),abs(sum(fprof_1(35:55,2:2:end))),'*'); % should use the abs values i think 
hold ;                                                         % i think the fsems has phase artifact
plot(TE(1:2:end),abs(sum(fprof_1(35:55,1:2:end))),'r*'); % possibly only fit even echoes


% Now do the fitting:
decay = sum(real(fprof_1(35:55,:)));
decay = decay(2:2:end);
TE = TE(2:2:end);

Mo_guess = 5*decay(1);
T2_guess = 0.050;

LB = [100, 0.0005];
UB = [10e8, 0.3];

guess0 = [Mo_guess; T2_guess];
guess = lsqnonlin('ps_lsq', ...
    guess0, LB, UB, ...
    [],...
    TE, ...
    decay);

M0 = guess(1)
T2 = guess(2)

m = ps_lsq(guess,TE);
hold on
plot(TE,m,'r')
hold off