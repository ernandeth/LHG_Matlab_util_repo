% creating ASL timing schedule that could be used to fit transit time
Nframes = 50

tmp = linspace(0.08, 3, Nframes/2);
t_tag = zeros(Nframes,1 );
t_tag(1:2:end) = tmp;
t_tag(2:2:end) = tmp;

t_delay =  ones(size(t_tag)) * 0.08;
isLabel = ones(size(t_tag));
isLabel(1:2:end) = 0;
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = 5 - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

!mkdir bat_schedule
save bat_schedule/t_tags.txt t_tag -ascii
save bat_schedule/t_adjusts.txt t_adjust -ascii
save bat_schedule/t_delays.txt t_delay -ascii
save bat_schedule/labelcontrol.txt isLabel -ascii


dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

%[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);
[dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);

subplot(311)
plot(dict');

d2 = dict(:, 2:2:end) - dict(:, 1:2:end); plot(d2');

dict = dict(:,10:end);

subplot(312)
plot(dict');

xc = corrcoef(abs(dict'));
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

G = mean(abs(xc(xc>0)))

subplot(313)
imagesc(xc)
caxis([0.8 1])
colorbar



%%  Now make a schedule for a saturation recovery experiment

Nframes = 10


t_tag = ones(6*Nframes,1 ) * 0.08;
t_delay =  ones(size(t_tag)) * 0.08;
isLabel = -ones(size(t_tag));
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = linspace(0.2, 5, Nframes)' ;

tmp=repmat(t_adjust,1,6)' ;
t_adjust = tmp(:);
t_adjust = t_adjust - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

!mkdir SR_schedule
save SR_schedule/t_tags.txt t_tag -ascii
save SR_schedule/t_adjusts.txt t_adjust -ascii
save SR_schedule/t_delays.txt t_delay -ascii
save SR_schedule/labelcontrol.txt isLabel -ascii


%

dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

subplot(311)
plot(dict');

d2 = dict(:, 2:2:end) - dict(:, 1:2:end); plot(d2');

dict = dict(:,10:end);

subplot(312)
plot(dict');

xc = corrcoef(abs(dict'));
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

G = mean(abs(xc(xc>0)))

subplot(313)
imagesc(xc)
caxis([0.8 1])
colorbar



%%  Now make a schedule for the perfusion part of the ASL experiment

Nframes = 50


t_tag = ones(Nframes,1 ) * 1.8;
t_delay =  ones(size(t_tag)) * 1.8;
isLabel = ones(size(t_tag));
isLabel(1:2:end) = 0;
isLabel(1:8) = 0;
t_aq = ones(size(t_tag)) * 0.035;
t_adjust = ones(size(t_tag)) * 4  ;
t_adjust = t_adjust - t_tag - t_delay - t_aq;

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

!mkdir CASL_schedule
save CASL_schedule/t_tags.txt t_tag -ascii
save CASL_schedule/t_adjusts.txt t_adjust -ascii
save CASL_schedule/t_delays.txt t_delay -ascii
save CASL_schedule/labelcontrol.txt isLabel -ascii


dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;

[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);

subplot(311)
plot(dict');

d2 = dict(:, 2:2:end) - dict(:, 1:2:end); plot(d2');

dict = dict(:,10:end);

subplot(312)
plot(dict');

xc = corrcoef(abs(dict'));
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

G = mean(abs(xc(xc>0)))

subplot(313)
imagesc(xc)
caxis([0.8 1])
colorbar


%% Make a schedule based on sine fluctuations.
Nframes = 150


t_tag = 2.0 * abs(sin(linspace(0,3*pi,Nframes)));
t_delay =  1.6 * abs(sin(linspace(0,5*pi,Nframes)));
isLabel = round(rand(size(t_delay)));

t_aq = ones(size(t_tag)) * 0.035;
t_adjust = ones(size(t_tag)) * 4  ;
t_adjust = t_adjust - t_tag - t_delay - t_aq;


plot(t_tag)
hold on
plot(t_delay)
stem(isLabel)
figure

timing_parms.Nlabel_group = 1;
timing_parms.t_tag = t_tag;
timing_parms.t_delay = t_delay;
timing_parms.t_adjust = t_adjust;
timing_parms.isLabel = isLabel;
timing_parms.order = 1;
timing_parms.t_aq = t_aq;

!mkdir sine_schedule
save sine_schedule/t_tags.txt t_tag -ascii
save sine_schedule/t_adjusts.txt t_adjust -ascii
save sine_schedule/t_delays.txt t_delay -ascii
save sine_schedule/labelcontrol.txt isLabel -ascii


%%
Npoints = 300;
timing_parms.t_tag = 1*abs(sinc(linspace(2, 0, Npoints))) + 0.05 ;
timing_parms.t_delay =  0.05 * ones(1,Npoints) ;
timing_parms.t_adjust = 0.05 * ones(1,Npoints) ;
timing_parms.isLabel = round(rand(1,Npoints))  ;
supercontrols = randperm(Npoints);
supercontrols=supercontrols(1:round(Npoints/4));
timing_parms.isLabel(supercontrols) = -1;
timing_parms.order = 1;
timing_parms.Nlabel_group = 1;
timing_parms.t_aq = ones(1,Npoints) * 0.035;

t_tag = timing_parms.t_tag;
t_delay = timing_parms.t_delay;
t_adjust = timing_parms.t_adjust;
isLabel = timing_parms.isLabel;

!mkdir sinc_schedule
save sinc_schedule/t_tags.txt t_tag -ascii
save sinc_schedule/t_adjusts.txt t_adjust -ascii
save sinc_schedule/t_delays.txt t_delay -ascii
save sinc_schedule/labelcontrol.txt isLabel -ascii

%%

dict_phys_parms.r1tis = 1./linspace(0.5, 3, 3);;
dict_phys_parms.flip =   deg2rad([50, 60, 70]);
dict_phys_parms.bat =  linspace(0.5, 4, 3);
dict_phys_parms.f =  linspace(0,100,3) / 6000;
dict_phys_parms.cbva = [linspace(0, 0.03, 3)];% , linspace(0.035,1,5)];
dict_phys_parms.kfor = 0.02;
dict_phys_parms.Disp = 20;
dict_phys_parms.mtis0 = 1;





%[dict, parms] = gen_flex_dictionary_150521 (timing_parms, dict_phys_parms);
[dict, parms] = gen_flex_dictionary_160525 (timing_parms, dict_phys_parms);

subplot(311)
plot(dict');

d2 = dict(:, 2:2:end) - dict(:, 1:2:end); plot(d2');

%dict = dict(:,10:end);

subplot(312)
plot(dict');

xc = corrcoef(abs(dict'));
% I only care about half of the matrix
for n=1:size(xc,1)
    xc(n, 1:n) = 0;
end

G = mean(abs(xc(xc>0)))

subplot(313)
imagesc(xc)
caxis([0.8 1])
colorbar


