
delta_guess = 1;
Tau_guess = 3;
amp_guess = 20;

guess0 = [  delta_guess Tau_guess amp_guess ];
guess_max = [3   6   50];
guess_min = [-1    1    0.0001];

optvar=optimset('lsqnonlin');
optvar.TolFun=1e-15;
optvar.TolX=1e-15;
optvar.MaxIter=600;
optvar.Display='iter';

t = linspace(0,99);

m = hrf_lsq( guess0 , t) + rand(1,100);


guess = lsqnonlin('hrf_lsq',...
    guess0, guess_min, guess_max, ...
    optvar,...
    t, ...
    m);


m2 = hrf_lsq( guess, t);
            
plot(t,m,'*', t,m2);
