subjects = ['020831js';
	'020831ld';
	'020915as';
	'020919rb';
	'020929jb';
	'030131mh';
	'030312jg';
	'030314en';
	'030314ep';
	'030314ms'];

paradigm = 'visual.pos.left'

delta_guess = 0;
Tau_guess = 2;

delta2_guess = 10;
Tau2_guess = 4;
amp2_guess = 0.5;

guess0 = [  delta_guess Tau_guess  delta2_guess Tau2_guess  amp2_guess];
guess_max = [0.5   6     10 12 1];
guess_min = [-1  1     3  1 0.01 ];

optvar=optimset('lsqnonlin');
optvar.TolFun=1e-15;
optvar.TolX=1e-15;
optvar.MaxIter=600;
%optvar.Display='iter';

t = [0:30]';
all_ts =[];

for count = 1:size(subjects,1)
	filename = sprintf('/net/stanley/data/hernan/inhib/%s/avg.%s.dat', subjects(count,:),paradigm);
	if (2==exist(filename))
		fprintf('\n reading %s',filename);
		ts = load(filename);
		ts = ts(:,1);
		ts = (ts-ts(1) ) / max(abs(ts));
		all_ts = [ all_ts ts ];
		plot(ts)
	else
		fprintf('\nNo file : %s', filename)
	end

end

mean_ts = mean(all_ts,2);
mean_ts = mean_ts/(max(abs(mean_ts)));

guess = lsqnonlin('hrf_lsq',...
    guess0, guess_min, guess_max, ...
    optvar,...
    t, ...
    mean_ts)


m2 = hrf_lsq( guess, t);
            
plot(t,mean_ts,'*', t,m2);
tmp = abs(m2);

peaktime = t(find(tmp == max(tmp)))

save timeseries.mat
