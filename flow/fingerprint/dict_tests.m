
tmp =       1-cos(linspace(0, 4*pi, 40));
% tmp = rand(1,40) *2; 
label_dur = [tmp tmp]; 
label_dur(1:2:end) = tmp;
label_dur(2:2:end) = tmp;

tmp = rand(1,40) *2; 
tmp =       1-cos(linspace(0, 4*pi, 40)) ;
PID = [tmp tmp]; 
PID(1:2:end) = tmp;
PID(2:2:end) = tmp;

tmp = 0.5* rand(1,40) + 0.7; 
tmp = 1-0.5*cos(linspace(0, 4*pi, 40)) ;
TR = [tmp tmp]; 
TR = TR + label_dur + PID ;

TR(1:2:end) = tmp;
TR(2:2:end) = tmp;

timing_parms.PID = PID';
timing_parms.label_dur = label_dur';
timing_parms.TR = TR';


phys.f =         0.01;
phys.f =         60 /6000;
phys.mtis0 =     1 ;
phys.cbva =      0.02 ;
phys.transit =   1.2 ;
phys.kfor =      0.2; % 1e-2 ;
phys.r1tis =     1/1.4  ;
phys.beta =      90*pi/180 ; % flip angle in radians
phys.Disp =      30;

dofigs = 1;
doSub=0;

obs = gen_signals_140618(phys, timing_parms, dofigs,doSub);
[dict, parms] = gen_dictionary_140618 (timing_parms);


% subtracted dictionary
sdict = dict(:,4:2:end) - dict(:, 3:2:end);

g = corr(sdict');
min(g(:))

subplot(311)
plot(sdict'); 

sdict = diff(dict, [],2);
sdict = sdict(:, 2:end);
g = corr(sdict');
min(g(:))
subplot(312)
plot(sdict'); 

dict2 = dict(:,2:end);
g = corr(dict2');
min(g(:))
subplot(313)
plot(dict2'); 


% noting that if flow is the only variable of interest, the best way to
% decorrelate thte dictionary entries is to pairwise subtract them
% and then remove the first point in the time series.
