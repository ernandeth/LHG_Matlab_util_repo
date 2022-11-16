function soln=bold_fit(data, SEC)
% points per second:
SEC=10;
t=0:max(size(data))-1;


% set initial guess and boundaries
Tau_guess=4*SEC;
delta_guess=2*SEC;
Amp_guess = 2;

guess0 = [   delta_guess Tau_guess Amp_guess ];
guess_max = [5*SEC   10*SEC   50];
guess_min = [0*SEC    0*SEC    -50];


% select options for search
optvar=optimset('lsqnonlin');
optvar.TolFun=1e-15;
optvar.TolX=1e-15;
optvar.MaxIter=600;
optvar.Display='off';   %'iter';


% generate some fake data
% data = BOLD_lsq( guess0 , t);
% data = data + 0.5*randn(size(data));

guess = lsqnonlin('BOLD_lsq',...
    guess0, guess_min, guess_max, ...
    optvar,...
    t, ...
    data);


data2= BOLD_lsq( guess, t);

plot(t/SEC,data)
hold on
plot(t/SEC,data2,'r')

soln=[ guess(1)/SEC guess(2)/SEC guess(3)];

return